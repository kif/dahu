#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* SubtractBuffer: Search for the equivalence of buffers, average them and subtract from sample signal.  
* SaxsAnalysis: Performs Guinier + Kratky + IFT, generates plots  
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "25/11/2024"
__status__ = "development"
__version__ = "0.2.1"

import os
import json
from math import log, pi
from collections import namedtuple
from urllib3.util import parse_url
from dahu.plugin import Plugin
from dahu.utils import fully_qualified_name
import logging
logger = logging.getLogger("bm29.subtract")
import numpy
try:
    import numexpr
except ImportError:
    logger.error("Numexpr is not installed, falling back on numpy's implementations")
    numexpr = None
import h5py
import pyFAI, pyFAI.azimuthalIntegrator
from pyFAI.containers import Integrate1dResult
from pyFAI.method_registry import IntegrationMethod
import freesas, freesas.cormap, freesas.invariants
from freesas.autorg import auto_gpa, autoRg, auto_guinier
from freesas.bift import BIFT
from scipy.optimize import minimize
from .common import Sample, Ispyb, get_equivalent_frames, cmp_float, get_integrator, KeyCache, \
                    polarization_factor, method, Nexus, get_isotime, SAXS_STYLE, NORMAL_STYLE, \
                    Sample, create_nexus_sample
from .ispyb import IspybConnector, NumpyEncoder

try:
    import memcache
except (ImportError, ModuleNotFoundError):
    memcache = None

NexusJuice = namedtuple("NexusJuice", "filename h5path npt unit q I sigma poni mask energy polarization method signal2d error2d normalization sample")


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
        self.output["memcached"] = self.send_to_memcached()
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
        self.sample_juice = self.read_nexus(self.sample_file)
        self.to_pyarch["basename"] = os.path.splitext(os.path.basename(self.sample_file))[0]
        try:
            self.create_nexus()
        except Exception as err:
            # try to register in test-mode
            if self.input.get("test_mode", True):
                try:
                    self.send_to_ispyb()
                except Exception as err2:
                    import traceback
                    self.log_warning(f"Processing failed and unable to send remaining data to ISPyB: {type(err2)} {err2}\n{traceback.format_exc(limit=10)}")
                raise(err)
        else:
            self.send_to_ispyb()        
        


    def validate_buffer(self, buffer_file):
        "Validate if a buffer is consitent with the sample, return some buffer_juice or None when unconsistent"
        buffer_juice = self.read_nexus(buffer_file)
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

    def create_nexus(self):
        nxs = self.nxs = Nexus(self.output_file, mode="w")
        entry_grp = nxs.new_entry("entry", self.input.get("plugin_name", "dahu"),
                                  title='BioSaxs buffer subtraction',
                                  force_time=get_isotime())
        nxs.h5.attrs["default"] = entry_grp.name

    # Configuration
        cfg_grp = nxs.new_class(entry_grp, "configuration", "NXnote")
        cfg_grp.create_dataset("data", data=json.dumps(self.input, indent=2, separators=(",\r\n", ":\t")))
        cfg_grp.create_dataset("format", data="text/json")

    # Process 0: Measurement group
        input_grp = nxs.new_class(entry_grp, "0_measurement", "NXcollection")
        input_grp["sequence_index"] = 0
        rel_path = os.path.relpath(os.path.abspath(self.sample_file), os.path.dirname(os.path.abspath(self.output_file)))
        input_grp["sample"] = h5py.ExternalLink(rel_path, self.sample_juice.h5path)

        for idx, buffer_file in enumerate(self.buffer_files):
            buffer_juice = self.validate_buffer(buffer_file)
            if buffer_juice is not None:
                rel_path = os.path.relpath(os.path.abspath(buffer_file), os.path.dirname(os.path.abspath(self.output_file)))
                input_grp["buffer_%i" % idx] = h5py.ExternalLink(rel_path, buffer_juice.h5path)
                self.buffer_juices.append(buffer_juice)

        # Sample: outsourced !
        create_nexus_sample(nxs, entry_grp, self.sample_juice.sample)

    # Process 1: CorMap
        cormap_grp = nxs.new_class(entry_grp, "1_correlation_mapping", "NXprocess")
        cormap_grp["sequence_index"] = 1
        cormap_grp["program"] = "freesas.cormap"
        cormap_grp["version"] = freesas.version
        cormap_grp["date"] = get_isotime()
        cormap_data = nxs.new_class(cormap_grp, "results", "NXdata")
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
        cormap_grp.attrs["default"] = cormap_data.name

    # Process 2: Image processing: subtraction with standard deviation
        average_grp = nxs.new_class(entry_grp, "2_buffer_subtraction", "NXprocess")
        average_grp["sequence_index"] = 2
        average_grp["program"] = fully_qualified_name(self.__class__)
        average_grp["version"] = __version__
        average_data = nxs.new_class(average_grp, "results", "NXdata")
        average_data.attrs["SILX_style"] = SAXS_STYLE
        average_data.attrs["signal"] = "intensity_normed"
    # Stage 2 processing

        # Nota: This formula takes into account the number of input frames in each averaged buffer !        
        #  avg = Σdata / Σnorm
        #  var = sigma² = ΣV / Σnorm
        # TODO implement those math using numexpr:
        sum_signal = sum(self.buffer_juices[i].normalization * self.buffer_juices[i].signal2d for i in to_merge_idx)
        sum_variance = sum((self.buffer_juices[i].normalization * self.buffer_juices[i].error2d) ** 2 for i in to_merge_idx)
        sum_norm = sum(self.buffer_juices[i].normalization for i in to_merge_idx)
        sum_norm2 = sum(self.buffer_juices[i].normalization**2 for i in to_merge_idx)
        buffer_average = sum_signal / sum_norm
        buffer_variance = sum_variance / sum_norm2
        sample_average = self.sample_juice.signal2d
        sample_variance = self.sample_juice.error2d ** 2
        sub_average = sample_average - buffer_average
        sub_variance = sample_variance + buffer_variance
        sub_std = numpy.sqrt(sub_variance)

        int_avg_ds = average_data.create_dataset("intensity_normed",
                                                 data=numpy.ascontiguousarray(sub_average, dtype=numpy.float32),
                                                 **cmp_float)
        int_avg_ds.attrs["interpretation"] = "image"
        int_avg_ds.attrs["formula"] = "sample_signal - weighted_mean(buffer_signal_i)"
        int_std_ds = average_data.create_dataset("intensity_std",
                                                 data=numpy.ascontiguousarray(sub_std, dtype=numpy.float32),
                                                 **cmp_float)
        int_std_ds.attrs["interpretation"] = "image"
        int_std_ds.attrs["formula"] = "sqrt( sample_variance + weighted_mean(buffer_variance_i) )"
        int_std_ds.attrs["method"] = "quadratic sum of sample error and buffer errors"
        average_grp.attrs["default"] = average_data.name

        key_cache = KeyCache(self.sample_juice.npt, self.sample_juice.unit, self.sample_juice.poni, self.sample_juice.mask, self.sample_juice.energy)
        ai = get_integrator(key_cache)

        if self.ispyb.url and parse_url(self.ispyb.url).host:
            # we need to provide the sample record and the best_buffer so let's generate it
            # This is a waist of time & resources ...
            res1 = ai._integrate1d_ng(sample_average, key_cache.npt,
                                      variance=sample_variance,
                                      polarization_factor=self.sample_juice.polarization,
                                      unit=key_cache.unit,
                                      safe=False,
                                      method=self.sample_juice.method)
            self.to_pyarch["sample"] = res1
            # Integrate also the buffer
            res2 = ai._integrate1d_ng(buffer_average, key_cache.npt,
                                      variance=buffer_variance,
                                      polarization_factor=self.sample_juice.polarization,
                                      unit=key_cache.unit,
                                      safe=False,
                                      method=self.sample_juice.method)
            self.to_pyarch["buffer"] = res2

    # Process 3: Azimuthal integration of the subtracted image
        ai2_grp = nxs.new_class(entry_grp, "3_azimuthal_integration", "NXprocess")
        ai2_grp["sequence_index"] = 3
        ai2_grp["program"] = "pyFAI"
        ai2_grp["version"] = pyFAI.version
        ai2_grp["date"] = get_isotime()
        radial_unit, unit_name = str(key_cache.unit).split("_", 1)
        ai2_data = nxs.new_class(ai2_grp, "results", "NXdata")
        ai2_data.attrs["SILX_style"] = SAXS_STYLE
        ai2_data.attrs["title"] = "%s, subtracted" % self.sample_juice.sample.name
        ai2_data.attrs["signal"] = "I"
        ai2_data.attrs["axes"] = radial_unit
        ai2_grp.attrs["default"] = ai2_data.name
        cfg_grp = nxs.new_class(ai2_grp, "configuration", "NXnote")
        cfg_grp.create_dataset("data", data=json.dumps(ai.get_config(), indent=2, separators=(",\r\n", ": ")))
        cfg_grp.create_dataset("format", data="text/json")
        if os.path.exists(key_cache.poni):
            cfg_grp.create_dataset("file_name", data=key_cache.poni)
        pol_ds = cfg_grp.create_dataset("polarization_factor", data=polarization_factor)
        pol_ds.attrs["comment"] = "Between -1 and +1, 0 for circular"
        cfg_grp.create_dataset("integration_method", data=json.dumps(method.method._asdict()))

    # Stage 3 processing: azimuthal integration
        res3 = ai._integrate1d_ng(sub_average, key_cache.npt,
                                  variance=sub_variance,
                                  polarization_factor=self.sample_juice.polarization,
                                  unit=key_cache.unit,
                                  safe=False,
                                  method=self.sample_juice.method)

        ai2_q_ds = ai2_data.create_dataset(radial_unit,
                                           data=numpy.ascontiguousarray(res3.radial, dtype=numpy.float32))
        ai2_q_ds.attrs["units"] = unit_name
        radius_unit = "nm" if "nm" in unit_name else "Å"
        ai2_q_ds.attrs["long_name"] = f"Scattering vector q ({radius_unit}⁻¹)"

        ai2_int_ds = ai2_data.create_dataset("I", data=numpy.ascontiguousarray(res3.intensity, dtype=numpy.float32))
        ai2_std_ds = ai2_data.create_dataset("errors",
                                             data=numpy.ascontiguousarray(res3.sigma, dtype=numpy.float32))

        ai2_int_ds.attrs["interpretation"] = "spectrum"
        ai2_int_ds.attrs["units"] = "arbitrary"
        ai2_int_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"
        #  ai2_int_ds.attrs["uncertainties"] = "errors" #this does not work
        ai2_std_ds.attrs["interpretation"] = "spectrum"
        ai2_int_ds.attrs["units"] = "arbitrary"

        if self.ispyb.url and parse_url(self.ispyb.url).host:
            self.to_pyarch["subtracted"] = res3
        # Export this to the output JSON
        #self.output["q"] = res3.radial
        #self.output["I"] = res3.intensity
        #self.output["std"] = res3.sigma
        self.to_memcached["q"] = res3.radial
        self.to_memcached["I"] = res3.intensity
        self.to_memcached["std"] = res3.sigma

        #  Finally declare the default entry and default dataset ...
        #  overlay the BIFT fitted data on top of the scattering curve
        entry_grp.attrs["default"] = ai2_data.name

    # Process 4: Guinier analysis
        guinier_grp = nxs.new_class(entry_grp, "4_Guinier_analysis", "NXprocess")
        guinier_grp["sequence_index"] = 4
        guinier_grp["program"] = "freesas.autorg"
        guinier_grp["version"] = freesas.version
        guinier_grp["date"] = get_isotime()
        guinier_autorg = nxs.new_class(guinier_grp, "autorg", "NXcollection")
        guinier_gpa = nxs.new_class(guinier_grp, "gpa", "NXcollection")
        guinier_guinier = nxs.new_class(guinier_grp, "guinier", "NXcollection")
        guinier_data = nxs.new_class(guinier_grp, "results", "NXdata")
        guinier_data.attrs["SILX_style"] = NORMAL_STYLE
        guinier_data.attrs["title"] = "Guinier analysis"
    # Stage4 processing: autorg and auto_gpa
        sasm = numpy.vstack((res3.radial, res3.intensity, res3.sigma)).T

        try:
            gpa = auto_gpa(sasm)
        except Exception as error:
            guinier_gpa["Failed"] = "%s: %s" % (error.__class__.__name__, error)
            gpa = None
        else:
            #  "Rg sigma_Rg I0 sigma_I0 start_point end_point quality aggregated"
            guinier_gpa["Rg"] = gpa.Rg
            guinier_gpa["Rg"].attrs["unit"] = radius_unit
            guinier_gpa["Rg_error"] = gpa.sigma_Rg
            guinier_gpa["Rg_error"].attrs["unit"] = radius_unit
            guinier_gpa["I0"] = gpa.I0
            guinier_gpa["I0_error"] = gpa.sigma_I0
            guinier_gpa["start_point"] = gpa.start_point
            guinier_gpa["end_point"] = gpa.end_point
            guinier_gpa["quality"] = gpa.quality
            guinier_gpa["aggregated"] = gpa.aggregated

        try:
            guinier = auto_guinier(sasm)
        except Exception as error:
            guinier_guinier["Failed"] = "%s: %s" % (error.__class__.__name__, error)
            guinier = None
        else:
            #  "Rg sigma_Rg I0 sigma_I0 start_point end_point quality aggregated"
            guinier_guinier["Rg"] = guinier.Rg
            guinier_guinier["Rg"].attrs["unit"] = radius_unit
            guinier_guinier["Rg_error"] = guinier.sigma_Rg
            guinier_guinier["Rg_error"].attrs["unit"] = radius_unit
            guinier_guinier["I0"] = guinier.I0
            guinier_guinier["I0_error"] = guinier.sigma_I0
            guinier_guinier["start_point"] = guinier.start_point
            guinier_guinier["end_point"] = guinier.end_point
            guinier_guinier["quality"] = guinier.quality
            guinier_guinier["aggregated"] = guinier.aggregated
            guinier_guinier["qₘᵢₙ·Rg"] = guinier.Rg * res3.radial[guinier.start_point]
            guinier_guinier["qₘₐₓ·Rg"] = guinier.Rg * res3.radial[guinier.end_point - 1]

        try:
            autorg = autoRg(sasm)
        except Exception as err:
            guinier_autorg["Failed"] = "%s: %s" % (err.__class__.__name__, err)
            autorg = None
        else:
            if autorg.Rg < 0:
                guinier_autorg["Failed"] = "No Guinier region found with this algorithm"
                autorg = None
            else:
                # "Rg sigma_Rg I0 sigma_I0 start_point end_point quality aggregated"
                guinier_autorg["Rg"] = autorg.Rg
                guinier_autorg["Rg"].attrs["unit"] = radius_unit
                guinier_autorg["Rg_error"] = autorg.sigma_Rg
                guinier_autorg["Rg_error"].attrs["unit"] = radius_unit
                guinier_autorg["I0"] = autorg.I0
                guinier_autorg["I0_error"] = autorg.sigma_I0
                guinier_autorg["start_point"] = autorg.start_point
                guinier_autorg["end_point"] = autorg.end_point
                guinier_autorg["quality"] = autorg.quality
                guinier_autorg["aggregated"] = autorg.aggregated
                guinier_autorg["qₘᵢₙ·Rg"] = autorg.Rg * res3.radial[autorg.start_point]
                guinier_autorg["qₘₐₓ·Rg"] = autorg.Rg * res3.radial[autorg.end_point - 1]

        #  take one of the fits
        if guinier:
            guinier_data["source"] = "auto_guinier"
        elif autorg:
            guinier = autorg
            guinier_data["source"] = "autorg"
        elif gpa:
            guinier = gpa
            guinier_data["source"] = "gpa"
        else:
            guinier = None
            guinier_data["source"] = "None"

        if self.ispyb.url and parse_url(self.ispyb.url).host:
            self.to_pyarch["guinier"] = guinier

    # Stage #4 Guinier plot generation:

        q, I, err = sasm.T[:3]
        mask = (I > 0) & numpy.isfinite(I) & (q > 0) & numpy.isfinite(q)
        if err is not None:
            mask &= (err > 0.0) & numpy.isfinite(err)
        mask = mask.astype(bool)
        if guinier:
            self.Rg = guinier.Rg
            self.I0 = guinier.I0
#             first_point = guinier.start_point
#             last_point = guinier.end_point
            intercept = numpy.log(self.I0)
            slope = -self.Rg ** 2 / 3.0
            end = numpy.where(q > 1.5 / self.Rg)[0][0]
            mask[end:] = False

        q2 = q[mask] ** 2
        logI = numpy.log(I[mask])
        dlogI = err[mask] / logI
        q2_ds = guinier_data.create_dataset("q2", data=q2.astype(numpy.float32))
        q2_ds.attrs["unit"] = radius_unit + "⁻²"
        q2_ds.attrs["long_name"] = "q² (%s⁻²)" % radius_unit
        q2_ds.attrs["interpretation"] = "spectrum"
        lnI_ds = guinier_data.create_dataset("logI", data=logI.astype(numpy.float32))
        lnI_ds.attrs["long_name"] = "log(I)"
        lnI_ds.attrs["interpretation"] = "spectrum"
        erI_ds = guinier_data.create_dataset("errors", data=dlogI.astype(numpy.float32))
        erI_ds.attrs["interpretation"] = "spectrum"

        if guinier:
            guinier_data["fit"] = intercept + slope * q2
            guinier_data["fit"].attrs["slope"] = slope
            guinier_data["fit"].attrs["intercept"] = intercept

        guinier_data_attrs = guinier_data.attrs
        guinier_data_attrs["signal"] = "logI"
        guinier_data_attrs["axes"] = "q2"
        guinier_data_attrs["auxiliary_signals"] = "fit"
        guinier_grp.attrs["default"] = guinier_data.name
        if guinier is None:
            entry_grp.attrs["default"] = ai2_data.name
            self.log_error("No Guinier region found, data of dubious quality", do_raise=True)

    # Process 5: Kratky plot
        kratky_grp = nxs.new_class(entry_grp, "5_dimensionless_Kratky_plot", "NXprocess")
        kratky_grp["sequence_index"] = 5
        kratky_grp["program"] = "freesas.autorg"
        kratky_grp["version"] = freesas.version
        kratky_grp["date"] = get_isotime()
        kratky_data = nxs.new_class(kratky_grp, "results", "NXdata")
        kratky_data.attrs["SILX_style"] = NORMAL_STYLE
        kratky_data.attrs["title"] = "Dimensionless Kratky plots"
        kratky_grp.attrs["default"] = kratky_data.name

    # Stage #5 Kratky plot generation:
        Rg = guinier.Rg
        I0 = guinier.I0
        xdata = q * Rg
        ydata = xdata * xdata * I / I0
        dy = xdata * xdata * err / I0
        qRg_ds = kratky_data.create_dataset("qRg", data=xdata.astype(numpy.float32))
        qRg_ds.attrs["interpretation"] = "spectrum"
        qRg_ds.attrs["long_name"] = "q·Rg (unit-less)"
        
        #Nota the "/" hereafter is chr(8725), the division sign and not the usual slash
        k_ds = kratky_data.create_dataset("q2Rg2I∕I0", data=ydata.astype(numpy.float32)) 
        k_ds.attrs["interpretation"] = "spectrum"
        k_ds.attrs["long_name"] = "q²Rg²I(q)/I₀"
        ke_ds = kratky_data.create_dataset("errors", data=dy.astype(numpy.float32))
        ke_ds.attrs["interpretation"] = "spectrum"
        kratky_data_attrs = kratky_data.attrs
        kratky_data_attrs["signal"] = "q2Rg2I∕I0"
        kratky_data_attrs["axes"] = "qRg"

    # stage 6: Rambo-Tainer invariant
        rti_grp = nxs.new_class(entry_grp, "6_invariants", "NXprocess")
        rti_grp["sequence_index"] = 6
        rti_grp["program"] = "freesas.invariants"
        rti_grp["version"] = freesas.version
        rti_data = nxs.new_class(rti_grp, "results", "NXdata")
        # average_data.attrs["SILX_style"] = SAXS_STYLE
        # average_data.attrs["signal"] = "intensity_normed"
        # Rambo_Tainer
        rti = freesas.invariants.calc_Rambo_Tainer(sasm, guinier)
        Vc_ds = rti_data.create_dataset("Vc", data=rti.Vc)
        Vc_ds.attrs["unit"] = "nm²"
        Vc_ds.attrs["formula"] = "Rambo-Tainer: Vc = I₀/(sum_q qI(q) dq)"
        sigma_Vc_ds = rti_data.create_dataset("Vc_error", data=rti.sigma_Vc)
        sigma_Vc_ds.attrs["unit"] = "nm²"

        Qr_ds = rti_data.create_dataset("Qr", data=rti.Qr)
        Qr_ds.attrs["unit"] = "nm"
        Qr_ds.attrs["formula"] = "Rambo-Tainer: Qr = Vc/Rg"
        sigma_Qr_ds = rti_data.create_dataset("Qr_error", data=rti.sigma_Qr)
        sigma_Qr_ds.attrs["unit"] = "nm"

        mass_ds = rti_data.create_dataset("mass", data=rti.mass)
        mass_ds.attrs["unit"] = "kDa"
        mass_ds.attrs["formula"] = "Rambo-Tainer: mass = (Qr/ec)^(1/k)"
        sigma_mass_ds = rti_data.create_dataset("mass_error", data=rti.sigma_mass)
        sigma_mass_ds.attrs["unit"] = "kDa"

        self.Vc = rti.Vc
        self.mass = rti.mass

        volume = self.to_pyarch["volume"] = freesas.invariants.calc_Porod(sasm, guinier)
        volume_ds = rti_data.create_dataset("volume", data=volume)
        volume_ds.attrs["unit"] = "nm³"
        volume_ds.attrs["formula"] = "Porod: V = 2*π²I₀²/(sum_q I(q)q² dq)"
        self.to_pyarch["rti"] = rti

    # stage 7: Pair distribution function, what is the equivalent of datgnom
        bift_grp = nxs.new_class(entry_grp, "7_indirect_Fourier_transformation", "NXprocess")
        bift_grp["sequence_index"] = 6
        bift_grp["program"] = "freesas.bift"
        bift_grp["version"] = freesas.version
        bift_grp["date"] = get_isotime()
        bift_data = nxs.new_class(bift_grp, "results", "NXdata")
        bift_data.attrs["SILX_style"] = NORMAL_STYLE
        bift_data.attrs["title"] = "Pair distance distribution function p(r)"

        cfg_grp = nxs.new_class(bift_grp, "configuration", "NXcollection")
    # Process stage7, i.e. perform the IFT
        try:
            bo = BIFT(q, I, err)
            cfg_grp["Rg"] = guinier.Rg
            # Pretty limited quality as we have real time constrains
            cfg_grp["npt"] = npt = 64
            cfg_grp["Dmax÷Rg"] = 3
            Dmax = bo.set_Guinier(guinier, Dmax_over_Rg=3)
            # Pretty limited quality as we have real time constrains

            # First scan on alpha:
            cfg_grp["alpha_sup"] = alpha_max = bo.guess_alpha_max(npt)
            cfg_grp["alpha_inf"] = 1 / alpha_max
            cfg_grp["alpha_scan_steps"] = 11

            key = bo.grid_scan(Dmax, Dmax, 1,
                               1.0 / alpha_max, alpha_max, 11, npt)
            Dmax, alpha = key[:2]
            # Then scan on Dmax:
            cfg_grp["Dmax_sup"] = guinier.Rg * 4
            cfg_grp["Dmax_inf"] = guinier.Rg * 2
            cfg_grp["Dmax_scan_steps"] = 5
            key = bo.grid_scan(guinier.Rg * 2, guinier.Rg * 4, 5,
                               alpha, alpha, 1, npt)
            Dmax, alpha = key[:2]
            if bo.evidence_cache[key].converged:
                bo.update_wisdom()
                use_wisdom = True
            else:
                use_wisdom = False
            res = minimize(bo.opti_evidence, (Dmax, log(alpha)), args=(npt, use_wisdom), method="powell")
            cfg_grp["Powell_steps"] = res.nfev
            cfg_grp["Monte-Carlo_steps"] = 0
        except Exception as error:
            bift_grp["Failed"] = "%s: %s" % (error.__class__.__name__, error)
            bo = None
        else:
            stats = bo.calc_stats()
            bift_grp["alpha"] = stats.alpha_avg
            bift_grp["alpha_error"] = stats.alpha_std
            self.Dmax = bift_grp["Dmax"] = stats.Dmax_avg
            bift_grp["Dmax_error"] = stats.Dmax_std
            bift_grp["S0"] = stats.regularization_avg
            bift_grp["S0_error"] = stats.regularization_std
            bift_grp["Chi2r"] = stats.chi2r_avg
            bift_grp["Chi2r_error"] = stats.chi2r_std
            bift_grp["logP"] = stats.evidence_avg
            bift_grp["logP_error"] = stats.evidence_std
            bift_grp["Rg"] = stats.Rg_avg
            bift_grp["Rg_error"] = stats.Rg_std
            bift_grp["I0"] = stats.I0_avg
            bift_grp["I0_error"] = stats.I0_std
            # Now the plot:
            r_ds = bift_data.create_dataset("r", data=stats.radius.astype(numpy.float32))
            r_ds.attrs["interpretation"] = "spectrum"

            r_ds.attrs["unit"] = radius_unit
            r_ds.attrs["long_name"] = "radius r(%s)" % radius_unit
            p_ds = bift_data.create_dataset("p(r)", data=stats.density_avg.astype(numpy.float32))
            p_ds.attrs["interpretation"] = "spectrum"
            bift_data["errors"] = stats.density_std
            bift_data.attrs["signal"] = "p(r)"
            bift_data.attrs["axes"] = "r"

            r = stats.radius
            T = numpy.outer(q, r / pi)
            T = (4 * pi * (r[-1] - r[0]) / (len(r) - 1)) * numpy.sinc(T)
            bift_ds = ai2_data.create_dataset("BIFT", data=T.dot(stats.density_avg).astype(numpy.float32))
            bift_ds.attrs["interpretation"] = "spectrum"
            ai2_data.attrs["auxiliary_signals"] = "BIFT"
            bift_grp.attrs["default"] = bift_data.name
            self.to_pyarch["bift"] = stats

    @staticmethod
    def read_nexus(filename):
        "return some NexusJuice from a HDF5 file "
        with Nexus(filename, "r") as nxsr:
            entry_grp = nxsr.get_entries()[0]
            h5path = entry_grp.name
            nxdata_grp = nxsr.h5[entry_grp.attrs["default"]]
            signal = nxdata_grp.attrs["signal"]
            axis = nxdata_grp.attrs["axes"]
            I = nxdata_grp[signal][()]
            q = nxdata_grp[axis][()]
            sigma = nxdata_grp["errors"][()]
            npt = len(q)
            unit = pyFAI.units.to_unit(axis + "_" + nxdata_grp[axis].attrs["units"])
            integration_grp = nxdata_grp.parent
            poni = str(integration_grp["configuration/file_name"][()]).strip()
            if not os.path.exists(poni):
                poni = str(integration_grp["configuration/data"][()]).strip()
            polarization = integration_grp["configuration/polarization_factor"][()]
            method = IntegrationMethod.select_method(**json.loads(integration_grp["configuration/integration_method"][()]))[0]
            instrument_grp = nxsr.get_class(entry_grp, class_type="NXinstrument")[0]
            detector_grp = nxsr.get_class(instrument_grp, class_type="NXdetector")[0]
            mask = detector_grp["pixel_mask"].attrs["filename"]
            mono_grp = nxsr.get_class(instrument_grp, class_type="NXmonochromator")[0]
            energy = mono_grp["energy"][()]
            img_grp = nxsr.get_class(entry_grp["3_time_average"], class_type="NXdata")[0]
            image2d = img_grp["intensity_normed"][()]
            error2d = img_grp["intensity_std"][()]
            norm =  img_grp["normalization"][()] if "normalization" in img_grp else None
            # Read the sample description:
            sample_grp = nxsr.get_class(entry_grp, class_type="NXsample")[0]
            sample_name = sample_grp.name

            buffer = sample_grp["buffer"][()] if "buffer" in sample_grp else ""
            concentration = sample_grp["concentration"][()] if "concentration" in sample_grp else ""
            description = sample_grp["description"][()] if "description" in sample_grp else ""
            hplc = sample_grp["hplc"][()] if "hplc" in sample_grp else ""
            temperature = sample_grp["temperature"][()] if "temperature" in sample_grp else ""
            temperature_env = sample_grp["temperature_env"][()] if "temperature_env" in sample_grp else ""
            sample = Sample(sample_name, description, buffer, concentration, hplc, temperature_env, temperature)

        return NexusJuice(filename, h5path, npt, unit, q, I, sigma, poni, mask, energy, polarization, method, image2d, error2d, norm, sample)

    def send_to_ispyb(self):
        if self.ispyb.url and parse_url(self.ispyb.url).host:
            ispyb = IspybConnector(*self.ispyb)
            ispyb.send_subtracted(self.to_pyarch)
            self.to_pyarch["experiment_type"]="sampleChanger"
            self.to_pyarch["sample"] = self.sample_juice.sample
            ispyb.send_icat(data=self.to_pyarch)
        else:
            self.log_warning("Not sending to ISPyB: no valid URL %s" % self.ispyb.url)


    def send_to_memcached(self):
        "Send the content of self.to_memcached to the storage"
        keys = {}
        rc = {}
        if memcache is not None:
            mc = memcache.Client([('stanza', 11211)])
            key_base = self.output_file
            for k in sorted(self.to_memcached.keys(), key=lambda i:self.to_memcached[i].nbytes):
                key = key_base + "_" + k
                keys[k] = key
                value = json.dumps(self.to_memcached[k], cls=NumpyEncoder)
                rc[k] = mc.set(key, value)
        self.log_warning(f"Return codes for memcached {rc}")
        return keys

