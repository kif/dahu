"""Data Analysis plugin for BM29: BioSaxs

* SubtractBuffer: Search for the equivalence of buffers, average them and subtract from sample signal.
* SaxsAnalysis: Performs Guinier + Kratky + IFT, generates plots
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "05/10/2026"
__status__ = "development"
__version__ = "0.4.0"

import copy
import json
import logging
import os
import posixpath
import zipfile
from typing import NamedTuple

import freesas
import freesas.cormap
import h5py
import numpy
import pyFAI
import pyFAI.integrator.azimuthal
from freesas.app.extract_ascii import write_ascii
from pyFAI.containers import ErrorModel, Integrate1dResult
from pyFAI.method_registry import IntegrationMethod
from urllib3.util import parse_url

from dahu.plugin import Plugin
from dahu.utils import fully_qualified_name

from .analysis import saxs_analysis
from .common import (
    NORMAL_STYLE,
    SAXS_STYLE,
    Ispyb,
    Sample,
    SequenceIndex,
    create_nexus_sample,
    get_equivalent_frames,
    str_,
)
from .icat import send_icat
from .ispyb import IspybConnector, NumpyEncoder
from .memcached import to_memcached
from .nexus import Nexus, get_isotime

logger = logging.getLogger("bm29.subtract")
try:
    import numexpr
except ImportError:
    logger.error("Numexpr is not installed, falling back on numpy's implementations")
    numexpr = None


class NexusJuice(NamedTuple):
    filename: str
    h5path: str
    npt: int
    unit: str
    q: numpy.ndarray
    I: numpy.ndarray
    sigma: numpy.ndarray
    poni:str
    mask: numpy.ndarray
    energy: float
    polarization: float
    method: tuple
    sum_signal: numpy.ndarray
    sum_variance: numpy.ndarray
    sum_normalization: numpy.ndarray
    sample:str
    I_all: numpy.ndarray
    sigma_all: numpy.ndarray
    sum_normalization2: numpy.ndarray = None
    count: numpy.ndarray = None
    error_model: str = None

    @classmethod
    def read(cls, filename):
        """Extract some NexusJuice from a HDF5 file, alternative constructor

        :param filename: name of the file
        :return: NexusJuice instance
        """
        with Nexus(filename, "r") as nxsr:
            entry_grp = nxsr.get_entries()[0]
            h5path = entry_grp.name
            nxdata_grp = entry_grp[entry_grp.attrs["default"]]
            signal = str_(nxdata_grp.attrs["signal"])
            axis = nxdata_grp.attrs["axes"]
            if not isinstance(axis, str):  # list of axes, the radial one is the last
                axis = str_(axis[-1])
            I_ary = nxdata_grp[signal][()]
            q = nxdata_grp[axis][()]
            sigma = nxdata_grp["errors"][()]
            npt = len(q)
            unit = pyFAI.units.to_unit(axis + "_" + str_(nxdata_grp[axis].attrs["units"]))
            # Configuration of the azimuthal integration: next to the averaged data or in the integration step
            average_grp = nxdata_grp.parent
            if "configuration" in average_grp:
                config_grp = average_grp["configuration"]
            else:
                config_grp = entry_grp["1_integration/configuration"]
            poni = str_(config_grp["file_name"][()]).strip()
            if not os.path.exists(poni):
                poni = str_(config_grp["data"][()]).strip()
            polarization = config_grp["polarization_factor"][()]
            method = IntegrationMethod.select_method(**json.loads(config_grp["integration_method"][()]))[0]
            # Accumulators of the frames merged in the average, unavailable in former files
            if "accumulators" in average_grp:
                accumulators_grp = average_grp["accumulators"]
                sum_signal = accumulators_grp["sum_signal"][()]
                sum_variance = accumulators_grp["sum_variance"][()]
                sum_normalization = accumulators_grp["sum_normalization"][()]
                sum_normalization2 = accumulators_grp["sum_normalization2"][()] if "sum_normalization2" in accumulators_grp else None
                count = accumulators_grp["count"][()] if "count" in accumulators_grp else None
                error_model = str_(accumulators_grp.attrs.get("error_model", "VARIANCE"))
            else:
                sum_signal = sum_variance = sum_normalization = sum_normalization2 = count = error_model = None
            instrument_grp = nxsr.get_class(entry_grp, class_type="NXinstrument")[0]
            detector_grp = nxsr.get_class(instrument_grp, class_type="NXdetector")[0]
            mask = detector_grp["pixel_mask"].attrs["filename"]
            mono_grp = nxsr.get_class(instrument_grp, class_type="NXmonochromator")[0]
            energy = mono_grp["energy"][()]
            # Read the sample description:
            sample_grp = nxsr.get_class(entry_grp, class_type="NXsample")[0]
            sample_name = posixpath.basename(sample_grp.name)

            buffer = str_(sample_grp["buffer"][()] if "buffer" in sample_grp else "")
            concentration = sample_grp["concentration"][()] if "concentration" in sample_grp else ""
            description = str_(sample_grp["description"][()]) if "description" in sample_grp else ""
            hplc = str_(sample_grp["hplc"][()]) if "hplc" in sample_grp else ""
            temperature = sample_grp["temperature"][()] if "temperature" in sample_grp else ""
            temperature_env = sample_grp["temperature_env"][()] if "temperature_env" in sample_grp else ""
            sample = Sample(sample_name, description, buffer, concentration, hplc, temperature_env, temperature)

            if "1_integration" in entry_grp:
                I_all = entry_grp["1_integration/result/I"][()]
                sigma_all = entry_grp["1_integration/result/errors"][()]
            else:
                I_all = []
                sigma_all = []

        return cls(filename=filename,
                    h5path=h5path,
                    npt=npt,
                    unit=unit,
                    q=q,
                    I=I_ary,
                    sigma=sigma,
                    poni=poni,
                    mask=mask,
                    energy=energy,
                    polarization=polarization,
                    method=method,
                    sum_signal=sum_signal,
                    sum_variance=sum_variance,
                    sum_normalization=sum_normalization,
                    sample=sample,
                    I_all=I_all,
                    sigma_all=sigma_all,
                    sum_normalization2=sum_normalization2,
                    count=count,
                    error_model=error_model)

    def to_result(self):
        """Rebuild a pyFAI Integrate1dResult from the accumulators, i.e. to average several of them with `union`

        :return: Integrate1dResult instance
        """
        missing = [k for k in ("sum_signal", "sum_variance", "sum_normalization", "sum_normalization2", "count")
                   if getattr(self, k) is None]
        if missing:
            raise ValueError(f"Unable to rebuild an Integrate1dResult from {self.filename}: missing {', '.join(missing)}")
        sum_signal = numpy.asarray(self.sum_signal, dtype=numpy.float64)
        sum_variance = numpy.asarray(self.sum_variance, dtype=numpy.float64)
        sum_normalization = numpy.asarray(self.sum_normalization, dtype=numpy.float64)
        sum_normalization2 = numpy.asarray(self.sum_normalization2, dtype=numpy.float64)
        npt = sum_signal.size
        result = Integrate1dResult(numpy.asarray(self.q, dtype=numpy.float64),
                                   numpy.zeros(npt, dtype=numpy.float64),
                                   numpy.zeros(npt, dtype=numpy.float64))
        result._set_sum_signal(sum_signal)
        result._set_sum_variance(sum_variance)
        result._set_sum_normalization(sum_normalization)
        result._set_sum_normalization2(sum_normalization2)
        result._set_count(numpy.asarray(self.count, dtype=numpy.float64))
        result._set_sem(numpy.zeros(npt, dtype=numpy.float64))
        result._set_std(numpy.zeros(npt, dtype=numpy.float64))
        result._set_unit(self.unit)
        result._set_polarization_factor(self.polarization)
        result._set_method(self.method)
        result._set_error_model(ErrorModel[self.error_model or "VARIANCE"])
        return result.__recalculate_means__()


def save_zip(filename, sample_juice, buffer_juices):
    """Save a stack of I into a zipfile with each frames in a dat-file.

    :param filename: name of the zip-file
    :param sample_juice:
    :param buffer_juices: list of buffer juice
    :return: nothing
    """
    destz_sample = "sample/"
    destz_buffer = "buffer_%1i/"
    common = {"q": sample_juice.q}
    if sample_juice.sample:
        sample = sample_juice.sample
        if sample.name:
            common["sample"] = sample.name
            destz_sample += sample.name
        else:
            destz_sample += "sample"

        if sample.buffer:
            common["buffer"] = sample.buffer
            destz_buffer += str_(sample.buffer)
        else:
            destz_buffer += "buffer"

        if sample.temperature_env:
            common["storage temperature"] = sample.temperature_env

        if sample.temperature:
            common["exposure temperature"] = sample.temperature

        if sample.concentration:
            common["concentration"] = sample.concentration
    destz_sample +=  "_%04i.dat"
    destz_buffer +=  "_%04i.dat"
    res = {}
    # sample
    for idx, (i, s) in enumerate(zip(sample_juice.I_all, sample_juice.sigma_all)):
        r = copy.copy(common)
        r["I"] = i
        r["std"] = s
        res[destz_sample % idx] = r
    # buffers
    for buffer_idx, buffer in enumerate(buffer_juices):
        for idx, (i, s) in enumerate(zip(buffer.I_all, buffer.sigma_all)):
            r = copy.copy(common)
            r["I"] = i
            r["std"] = s
            res[destz_buffer % (buffer_idx, idx)] = r

    with zipfile.ZipFile(filename, "w") as z:
        for name, frame in res.items():
            z.writestr(name, write_ascii(frame))


class SubtractBuffer(Plugin):
    """Search for the equivalence of buffers, average them and subtract from sample signal.

        Typical JSON file:
    {
      "buffer_files": ["buffer_001.h5", "buffer_002.h5"],
      "sample_file": "sample.h5",
      "output_file": "subtracted.h5",
      "fidelity": 0.001,
      "ispyb": {
        "url": "http://ispyb.esrf.fr:1234",
        "login": "mx1234",
        "passwd": "secret",
        "pyarch": "/data/pyarch/mx1234/sample"
        "measurement_id": -1,
        "collection_id": -1
       },
      "wait_for": [jobid_buffer1, jobid_buffer2, jobid_sample],
      "plugin_name": "bm29.subtractbuffer"
    }
    """

    def __init__(self):
        Plugin.__init__(self)
        self.buffer_files = []
        self.sample_file = None
        self.nxs = None
        self.output_file = None
        self.npt = None
        self.poni = None
        self.mask = None
        self.energy = None
        self.sample_juice = None
        self.buffer_juices = []
        self.Rg = self.I0 = self.Dmax = self.Vc = self.mass = None
        self.ispyb = None
        self.to_pyarch = {}
        self.to_memcached = {}  # data to be shared via memcached
        self.seq = SequenceIndex(0)

    def setup(self, kwargs=None):
        logger.debug("SubtractBuffer.setup")
        Plugin.setup(self, kwargs)

        wait_for = self.input.get("wait_for")
        if wait_for:
            for job_id in wait_for:
                self.wait_for(job_id)

        self.sample_file = self.input.get("sample_file")
        if self.sample_file is None:
            self.log_error("No sample file provided", do_raise=True)
        self.output_file = self.input.get("output_file")
        if self.output_file is None:
            lst = list(os.path.splitext(self.sample_file))
            lst.insert(1, "-sub")
            dirname, basename = os.path.split("".join(lst))
            dirname = os.path.dirname(dirname)
#            dirname = os.path.join(dirname, "processed")
            dirname = os.path.join(dirname, "subtract")
            self.output_file = os.path.join(dirname, basename)
            if not os.path.isdir(dirname):
                try:
                    os.makedirs(dirname)
                except Exception as err:
                    self.log_warning(f"Unable to create dir {dirname}. {type(err)}: {err}")

        self.buffer_files = [os.path.abspath(fn) for fn in self.input.get("buffer_files", [])
                             if os.path.exists(fn)]
        #Manage gallery here
        dirname = os.path.dirname(self.output_file)
        gallery = os.path.join(dirname, "gallery")
        if not os.path.isdir(gallery):
            try:
                os.makedirs(gallery)
            except Exception as err:
                self.log_warning(f"Unable to create dir {gallery}. {type(err)}: {err}")
        ispydict = self.input.get("ispyb", {})
        ispydict["gallery"] = gallery
        self.ispyb = Ispyb._fromdict(ispydict)

    def teardown(self):
        Plugin.teardown(self)
        logger.debug("SubtractBuffer.teardown")

        # export the output file location
        self.output["output_file"] = self.output_file
        self.output["Rg"] = self.Rg
        self.output["I0"] = self.I0
        self.output["Dmax"] = self.Dmax
        self.output["Vc"] = self.Vc
        self.output["mass"] = self.mass

        #teardown everything else:
        if self.nxs is not None:
            self.nxs.close()
        self.sample_juice = None
        self.buffer_juices = []
        self.ispyb = None
        self.to_pyarch = None
        self.to_memcached = None

    def process(self):
        Plugin.process(self)
        logger.debug("SubtractBuffer.process")
        self.sample_juice = NexusJuice.read(self.sample_file)
        self.to_pyarch["basename"] = os.path.splitext(os.path.basename(self.sample_file))[0]
        try:
            self.create_nexus()
        except Exception:
            # try to register in test-mode
            if self.input.get("test_mode", True):
                try:
                    self.send_to_ispyb()
                except Exception as err2:
                    import traceback
                    self.log_warning(f"Processing failed and unable to send remaining data to ISPyB: {type(err2)} {err2}\n{traceback.format_exc(limit=10)}")
                raise
        else:
            self.send_to_ispyb()
            self.send_to_icat()
        self.output["memcached"] = self.send_to_memcached()


    def validate_buffer(self, buffer_file):
        "Validate if a buffer is consitent with the sample, return some buffer_juice or None when unconsistent"
        buffer_juice = NexusJuice.read(buffer_file)
        if self.sample_juice.npt != buffer_juice.npt:
            self.log_warning(f"Sample {buffer_file} differs in number of points, discarding")
            return
        if abs(self.sample_juice.q - buffer_juice.q).max() > 1e-6:
            self.log_warning(f"Sample {buffer_file} differs in q-position, discarding")
            return
        if self.sample_juice.poni != buffer_juice.poni:
            self.log_warning(f"Sample {buffer_file} differs in poni-file, discarding")
            return
        if self.sample_juice.mask != buffer_juice.mask:
            self.log_warning(f"Sample {buffer_file} differs in mask-file, discarding")
            return
        if self.sample_juice.polarization != buffer_juice.polarization:
            self.log_warning(f"Sample {buffer_file} differs in polarization factor, discarding")
            return
        if self.sample_juice.sample.buffer != buffer_juice.sample.buffer:
            self.log_warning(f"Sample {buffer_file} differs in buffer descsription, discarding")
            return
        if buffer_juice.sample.concentration:
            self.log_warning(f"Sample {buffer_file} concentration not null, discarding")
            return
        if "buffers" not in self.to_pyarch:
            buffers = self.to_pyarch["buffers"] = []
        else:
            buffers = self.to_pyarch["buffers"]
        buffers.append(Integrate1dResult(buffer_juice.q, buffer_juice.I, buffer_juice.sigma))
        return buffer_juice

    def average_buffers(self, indices):
        """Average the buffers which are equivalent, according to CorMap.

        An Integrate1dResult is rebuilt from the accumulators of each buffer and they are merged with `union`,
        i.e. sum of signal, normalization and variance:
            I = Σ sum_signalᵢ / Σ sum_normalizationᵢ
            σ = √(Σ sum_varianceᵢ) / Σ sum_normalizationᵢ

        :param indices: indices of the buffers to merge in self.buffer_juices
        :return: Integrate1dResult with the averaged buffer
        """
        results = [self.buffer_juices[i].to_result() for i in indices]
        if not results:
            self.log_error("No valid buffer to average", do_raise=True)
        average = results[0]
        for other in results[1:]:
            average = average.union(other)
        return average

    def create_nexus(self):
        nxs = self.nxs = Nexus(self.output_file, mode="w")
        entry_grp = nxs.new_entry("entry", self.input.get("plugin_name", "dahu"),
                                  title='BioSaxs buffer subtraction',
                                  force_time=get_isotime())
        nxs.h5.attrs["default"] = entry_grp.name.strip("/")

    # Configuration
        cfg_grp = nxs.new_class(entry_grp, "configuration", "NXnote")
        cfg_grp.create_dataset("data", data=json.dumps(self.input, indent=2, separators=(",\r\n", ":\t")))
        cfg_grp.create_dataset("format", data="text/json")

    # Process 0: Measurement group
        seq = self.seq()
        input_grp = nxs.new_class(entry_grp, f"{seq}_measurement", "NXcollection")
        input_grp["sequence_index"] = seq
        rel_path = os.path.relpath(os.path.abspath(self.sample_file), os.path.dirname(os.path.abspath(self.output_file)))
        input_grp["sample"] = h5py.ExternalLink(rel_path, self.sample_juice.h5path)

        for idx, buffer_file in enumerate(self.buffer_files):
            buffer_juice = self.validate_buffer(buffer_file)
            if buffer_juice is not None:
                rel_path = os.path.relpath(os.path.abspath(buffer_file), os.path.dirname(os.path.abspath(self.output_file)))
                input_grp[f"buffer_{idx}"] = h5py.ExternalLink(rel_path, buffer_juice.h5path)
                self.buffer_juices.append(buffer_juice)

        # Sample: outsourced !
        create_nexus_sample(nxs, entry_grp, self.sample_juice.sample)

        #save input curves as zipfile: TODO Check that this is is working:
        save_zip(os.path.splitext(self.output_file)[0]+".zip",
                 self.sample_juice,
                 self.buffer_juices)

    # Process 1: CorMap
        seq = self.seq()
        cormap_grp = nxs.new_class(entry_grp, f"{seq}_correlation_mapping", "NXprocess")
        cormap_grp["sequence_index"] = seq
        cormap_grp["program"] = "freesas.cormap"
        cormap_grp["version"] = freesas.version
        cormap_grp["date"] = get_isotime()
        cormap_data = nxs.new_class(cormap_grp, "result", "NXdata")
        cormap_data.attrs["SILX_style"] = NORMAL_STYLE
        cfg_grp = nxs.new_class(cormap_grp, "configuration", "NXcollection")

    # Cormap processing
        nb_frames = len(self.buffer_juices)
        count = numpy.empty((nb_frames, nb_frames), dtype=numpy.uint16)
        proba = numpy.empty((nb_frames, nb_frames), dtype=numpy.float32)
        for i in range(nb_frames):
            proba[i, i] = 1.0
            count[i, i] = 0
            for j in range(i):
                res = freesas.cormap.gof(self.buffer_juices[i].I, self.buffer_juices[j].I)
                proba[i, j] = proba[j, i] = res.P
                count[i, j] = count[j, i] = res.c
        fidelity = self.input.get("fidelity", 0)
        cfg_grp["fidelity_abs"] = fidelity
        cfg_grp["fidelity_rel"] = fidelity
        tomerge = get_equivalent_frames(proba, fidelity, fidelity)

        cormap_data.attrs["signal"] = "probability"
        cormap_ds = cormap_data.create_dataset("probability", data=proba)
        cormap_ds.attrs["interpretation"] = "image"
        cormap_ds.attrs["long_name"] = "Probability to be the same"

        count_ds = cormap_data.create_dataset("count", data=count)
        count_ds.attrs["interpretation"] = "image"
        count_ds.attrs["long_name"] = "Longest sequence where curves do not cross each other"

        to_merge_idx = numpy.arange(*tomerge, dtype=numpy.uint16)
        to_merge_ds = cormap_data.create_dataset("to_merge", data=to_merge_idx)
        # self.log_warning(f"to_merge: {tomerge}")
        to_merge_ds.attrs["long_name"] = "Index of equivalent frames"
        cormap_grp.attrs["default"] = posixpath.relpath(cormap_data.name, cormap_grp.name)

    # Process 2: Average the equivalent buffers and subtract them from the sample, in 1D
        seq = self.seq()
        sub_grp = nxs.new_class(entry_grp, f"{seq}_buffer_subtraction", "NXprocess")
        sub_grp["sequence_index"] = seq
        sub_grp["program"] = fully_qualified_name(self.__class__)
        sub_grp["version"] = __version__
        sub_grp["date"] = get_isotime()
        radial_unit, unit_name = str(self.sample_juice.unit).split("_", 1)
        radius_unit = "nm" if "nm" in unit_name else "Å"

    # Stage 2 processing
        # Average of the equivalent buffers, weighted by their normalization (i.e. the number of frames merged)
        buffer_average = self.average_buffers(to_merge_idx)
        # From now on, accumulators are no more needed: subtraction of curves with quadratic sum of errors
        q = numpy.asarray(self.sample_juice.q, dtype=numpy.float64)
        sample_I = numpy.asarray(self.sample_juice.I, dtype=numpy.float64)
        sample_sigma = numpy.asarray(self.sample_juice.sigma, dtype=numpy.float64)
        sub_I = sample_I - buffer_average.intensity
        sub_sigma = numpy.sqrt(sample_sigma**2 + buffer_average.sigma**2)
        subtracted = Integrate1dResult(q, sub_I, sub_sigma)

        buffer_data = nxs.new_class(sub_grp, "buffer_average", "NXdata")
        buffer_data.attrs["SILX_style"] = SAXS_STYLE
        buffer_data.attrs["title"] = f"Average of buffers {', '.join(str(i) for i in to_merge_idx)}"
        buffer_data.attrs["signal"] = "I"
        buffer_data.attrs["axes"] = radial_unit
        for name, data in ((radial_unit, q), ("I", buffer_average.intensity), ("errors", buffer_average.sigma)):
            buffer_data.create_dataset(name, data=numpy.ascontiguousarray(data, dtype=numpy.float32))
        buffer_data[radial_unit].attrs["units"] = unit_name
        buffer_data["I"].attrs["formula"] = "Σᵢ sum_signalᵢ / Σᵢ sum_normalizationᵢ"
        buffer_data["errors"].attrs["formula"] = "√(Σᵢ sum_varianceᵢ) / Σᵢ sum_normalizationᵢ"

        ai2_data = nxs.new_class(sub_grp, "result", "NXdata")
        ai2_data.attrs["SILX_style"] = SAXS_STYLE
        ai2_data.attrs["title"] = f"{self.sample_juice.sample.name}, subtracted"
        ai2_data.attrs["signal"] = "I"
        ai2_data.attrs["axes"] = radial_unit
        sub_grp.attrs["default"] = posixpath.relpath(ai2_data.name, sub_grp.name)

        ai2_q_ds = ai2_data.create_dataset(radial_unit, data=numpy.ascontiguousarray(q, dtype=numpy.float32))
        ai2_q_ds.attrs["units"] = unit_name
        ai2_q_ds.attrs["long_name"] = f"Scattering vector q ({radius_unit}⁻¹)"
        ai2_int_ds = ai2_data.create_dataset("I", data=numpy.ascontiguousarray(sub_I, dtype=numpy.float32))
        ai2_int_ds.attrs["interpretation"] = "spectrum"
        ai2_int_ds.attrs["units"] = "arbitrary"
        ai2_int_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"
        ai2_int_ds.attrs["formula"] = "I_sample - weighted_mean(I_buffer_i)"
        ai2_std_ds = ai2_data.create_dataset("errors", data=numpy.ascontiguousarray(sub_sigma, dtype=numpy.float32))
        ai2_std_ds.attrs["interpretation"] = "spectrum"
        ai2_std_ds.attrs["formula"] = "sqrt(sigma_sample² + sigma_buffer²)"
        ai2_std_ds.attrs["method"] = "quadratic sum of sample error and buffer errors"

        if self.ispyb.url and parse_url(self.ispyb.url).host:
            self.to_pyarch["sample"] = Integrate1dResult(q, sample_I, sample_sigma)
            self.to_pyarch["buffer"] = buffer_average
            self.to_pyarch["subtracted"] = subtracted
        self.to_memcached["q"] = q
        self.to_memcached["I"] = sub_I
        self.to_memcached["std"] = sub_sigma

        #  Finally declare the default entry and default dataset ...
        #  overlay the BIFT fitted data on top of the scattering curve
        entry_grp.attrs["default"] = posixpath.relpath(ai2_data.name, entry_grp.name)

    # Process 4-7: Guinier analysis, Kratky plot, invariants and BIFT
        sasm = numpy.vstack((q, sub_I, sub_sigma)).T
        analysis = saxs_analysis(nxs, entry_grp, sasm, radius_unit, self.seq, curve_data=ai2_data)
        guinier = analysis.guinier
        if self.ispyb.url and parse_url(self.ispyb.url).host:
            self.to_pyarch["guinier"] = guinier
        if guinier is None:
            entry_grp.attrs["default"] = posixpath.relpath(ai2_data.name, entry_grp.name)
            self.log_error("No Guinier region found, data of dubious quality", do_raise=True)
        self.Rg = guinier.Rg
        self.I0 = guinier.I0
        self.Vc = analysis.rti.Vc
        self.mass = analysis.rti.mass
        self.to_pyarch["volume"] = analysis.volume
        self.to_pyarch["rti"] = analysis.rti
        if analysis.bift is not None:
            self.Dmax = analysis.bift.Dmax_avg
            self.to_pyarch["bift"] = analysis.bift

    def send_to_ispyb(self):
        if self.ispyb.url and parse_url(self.ispyb.url).host:
            ispyb = IspybConnector(*self.ispyb)
            ispyb.send_subtracted(self.to_pyarch)
        else:
            self.log_warning(f"Not sending to ISPyB: no valid URL {self.ispyb.url}")

    def send_to_icat(self):
        to_icat = copy.copy(self.to_pyarch)
        to_icat["experiment_type"] = "sample-changer"
        if self.sample_juice is None:
            self.log_warning("Sample_juice is None in send_to_icat. Not sending garbage")
            return
        to_icat["sample"] = self.sample_juice.sample
        metadata = {"scanType": "subtraction"}
        raw = [os.path.dirname(os.path.abspath(i)) for i in self.buffer_files]
        raw.append(os.path.dirname(os.path.abspath(self.sample_file)))

        if not (self.ispyb.url and parse_url(self.ispyb.url).host):
            self.log_warning("Not sending to iCat: ISPyB metadata not valid")
            return

        return send_icat(sample=self.sample_juice.sample,
                         raw=raw,
                         path=os.path.dirname(os.path.abspath(self.output_file)),
                         data=to_icat,
                         dataset="subtraction",
                         gallery=self.ispyb.gallery or os.path.join(os.path.dirname(os.path.abspath(self.output_file)), "gallery"),
                         metadata=metadata)

    def send_to_memcached(self):
        "Send the content of self.to_memcached to the storage"
        dico = {}
        key_base = self.output_file
        for k in sorted(self.to_memcached.keys(), key=lambda i:self.to_memcached[i].nbytes):
            key = key_base + "_" + k
            dico[key] = json.dumps(self.to_memcached[k], cls=NumpyEncoder)

        return to_memcached(dico)

