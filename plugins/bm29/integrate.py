#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* IntegrateMultiframe: perform the integration of many frames contained in a HDF5 file and average them

"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "05/05/2025"
__status__ = "development"
__version__ = "0.3.0"

import os
import time
import json
import logging
import copy
from collections import namedtuple
from urllib3.util import parse_url
from dahu.plugin import Plugin
from dahu.factory import register
from dahu.utils import fully_qualified_name

import numpy
import h5py
import pyFAI
import pyFAI.integrator.azimuthal
import freesas
import freesas.cormap

from .common import Sample, Ispyb, get_equivalent_frames, cmp_int, cmp_float, get_integrator, KeyCache, \
                    method, polarization_factor, Nexus, get_isotime, SAXS_STYLE, NORMAL_STYLE, \
                    create_nexus_sample
from .ispyb import IspybConnector, NumpyEncoder
from .icat import send_icat
from .memcached import to_memcached


logger = logging.getLogger("bm29.integrate")
try:
    import numexpr
except ImportError:
    logger.error("Numexpr is not installed, falling back on numpy's implementations")
    numexpr = None

IntegrationResult = namedtuple("IntegrationResult", "radial intensity sigma")
CormapResult = namedtuple("CormapResult", "probability count tomerge")
AverageResult = namedtuple("AverageResult", "average deviation normalization")


@register
class IntegrateMultiframe(Plugin):
    """perform the integration of many frames contained in a HDF5 file and average them

    Input parameters:
    :param poni_file: configuration of the geometry
    :param input_file: path for the HDF5 file

    Typical JSON file:
    {
      "input_file": "/tmp/file1.h5",
      "output_file": "/tmp/file1.h5", # optional
      "frame_ids": [101, 102],
      "timestamps": [1580985678.47, 1580985678.58],
      "monitor_values": [1, 1.1],
      "storage_ring_current": [199.6, 199.5],
      "exposure_time": 0.1,
      "normalisation_factor": 1.0,
      "poni_file": "/tmp/example.poni",
      "mask_file": "/tmp/mask.edf",
      "npt": 1000,
      "energy": 12.0, #keV
      "fidelity_abs": 1e-5,
      "fidelity_rel": 1e-3,
      "hplc_mode": 0,
      "timeout": 10,
      "plugin_name": "bm29.integratemultiframe",
      "sample": {
        "name": "bsa",
        "description": "protein description like Bovine Serum Albumin",
        "buffer": "description of buffer, pH, ...",
        "concentration": 0,
        "hplc": "column name and chromatography conditions",
        "temperature": 20,
        "temperature_env": 20},
      "ispyb": {
        "url": "http://ispyb.esrf.fr:1234",
        "login": "mx1234",
        "passwd": "secret",
        "pyarch": "/data/pyarch/mx1234/sample",
        "measurement_id": -1,
        "collection_id": -1
       }
    }
    """

    COPY_IMAGES = False

    def __init__(self):
        Plugin.__init__(self)
        self.sample = None
        self.ispyb = None
        self.input_file = None
        self._input_frames = None
        self._start_time = "2000-01-01T00:00:00Z"
        self._end_time = "2000-01-01T00:00:00Z"
        self.output_file = None
        self.nxs = None
        self.nb_frames = None
        self.ai = None
        self.npt = 1000
        self.timeout = 10
        self.unit = pyFAI.units.to_unit("q_nm^-1")
        # self.polarization_factor = 0.9 --> constant
        self.poni = self.mask = None
        self.energy = None
        # self.method = IntegrationMethod.select_method(1, "no", "csr", "opencl")[0] -> constant
        self.monitor_values = None
        self.normalization_factor = None
        self.scale_factor = None
        self.to_pyarch = {}  # contains all the stuff to be sent to Ispyb and pyarch
        self.to_memcached = {}  # data to be shared via memcached

    def setup(self, kwargs=None):
        logger.debug("IntegrateMultiframe.setup")
        Plugin.setup(self, kwargs)

        self.sample = Sample._fromdict(self.input.get("sample", {}))
        if not self.sample.name:
            self.sample = Sample("Unknown sample", *self.sample[1:])

        self.timeout = self.input.get("timeout", self.timeout)
        self.input_file = self.input.get("input_file")
        if self.input_file is not None:
            self.wait_file(self.input_file, timeout=self.timeout)
        else:
            self.log_error(f"No valid input file provided {self.input_file}")

        self.to_pyarch["basename"] = os.path.splitext(os.path.basename(self.input_file))[0]
        if self.input_file is None:
            self.log_error("No input file provided", do_raise=True)
        self.output_file = self.input.get("output_file")
        if self.output_file is None:
            lst = list(os.path.splitext(self.input_file))
            lst.insert(1, "-integrate")
            dirname, basename = os.path.split("".join(lst))
            dirname = dirname.replace("RAW_DATA", "PROCESSED_DATA")
            dirname = os.path.dirname(dirname)
            # dirname = os.path.join(dirname, "processed")
            dirname = os.path.join(dirname, "integrate")
            self.output_file = os.path.join(dirname, basename)
            if not os.path.isdir(dirname):
                try:
                    os.makedirs(dirname)
                except Exception as err:
                    self.log_warning(f"Unable to create dir {dirname}. {type(err)}: {err}")

            self.log_warning(f"No output file provided, using: {self.output_file}")
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

        self.nb_frames = len(self.input.get("frame_ids", []))
        self.npt = self.input.get("npt", self.npt)
        self.unit = pyFAI.units.to_unit(self.input.get("unit", self.unit))
        self.poni = self.input.get("poni_file")
        if self.poni is None:
            self.log_error("No poni-file provided! aborting", do_raise=True)
        self.mask = self.input.get("mask_file")
        self.energy = self.input.get("energy")
        if self.energy is None:
            self.log_error("No energy provided! aborting", do_raise=True)
        else:
            self.energy = numpy.float32(self.energy)  # It is important to fix the datatype of the energy
        self.monitor_values = numpy.array(self.input.get("monitor_values", 1), dtype=numpy.float64)
        self.normalization_factor = float(self.input.get("normalization_factor", 1))
        self.scale_factor = float(self.input.get("exposure_time", 1)) / self.normalization_factor

    def teardown(self):
        Plugin.teardown(self)
        logger.debug("IntegrateMultiframe.teardown")
        # export the output file location
        self.output["output_file"] = self.output_file
        if self.nxs is not None:
            self.nxs.close()
        if self.ai is not None:
            self.ai = None
        # clean cache
        if self._input_frames is not None:
            self._input_frames = None
        self.monitor_values = None
        self.to_pyarch = None
        self.ispyb = None

    @property
    def input_frames(self):
        "For performance reasons, all frames are read in one bloc and cached, this returns a 3D numpy array"
        if self._input_frames is None:
            try:
                with Nexus(self.input_file, "r", timeout=self.timeout) as nxs:
                    entry = nxs.get_entries()[0]
                    if "measurement" in entry:
                        measurement = entry["measurement"]
                    else:
                        self.log_error("No measurement in entry: %s of data_file: %s" % (entry, self.input_file))
                    self._input_frames = measurement["data"][...]
                    try:
                        self._start_time = entry["start_time"][()]
                        self._end_time = entry["end_time"][()]
                    except Exception as err:
                        self.log_error("Unable to read time %s: %s" % (type(err), str(err)), do_raise=False)
            except Exception as err:
                self.log_error("Unable to read images %s: %s" % (type(err), str(err)), do_raise=True)
        return self._input_frames

    def process(self):
        "does the integration of a set of frames"
        logger.debug("IntegrateMultiframe.process")
        self.ai = get_integrator(KeyCache(self.npt, self.unit, self.poni, self.mask, self.energy))
        self.create_nexus()
        self.output["memcached"] = self.send_to_memcached()
        self.send_to_ispyb()
        #self.output["icat"] = 
        self.send_to_icat()

    def wait_file(self, filename, timeout=None):
        """Wait for a file to appear on a filesystem

        :param filename: name of the file
        :param timeout: time-out in seconds

        Raises an exception and ends the processing in case of missing file!
	"""
        timeout = self.timeout if timeout is None else timeout
        end_time = time.perf_counter() + timeout
        dirname = os.path.dirname(filename)
        while not os.path.isdir(dirname):
            if time.perf_counter() > end_time:
                self.log_error(f"Filename {filename} did not appear in {timeout} seconds")
            time.sleep(0.1)
            os.stat(os.path.dirname(dirname))

        while not os.path.exists(filename):
            if time.perf_counter() > end_time:
                self.log_error(f"Filename {filename} did not appear in {timeout} seconds")
            time.sleep(0.1)
            os.stat(dirname)

        while not os.stat(filename).st_size:
            if time.perf_counter() > end_time:
                self.log_error(f"Filename {filename} did not appear in {timeout} seconds")
            time.sleep(0.1)
            os.stat(dirname)
        appear_time = time.perf_counter() - end_time + timeout
        if appear_time > 1.0:
            self.log_warning(f"Filename {filename} took {appear_time:.3f}s to appear on filesystem!")

    def create_nexus(self):
        "create the nexus result file with basic structure"
        if not os.path.isdir(os.path.dirname(self.output_file)):
            os.makedirs(os.path.dirname(self.output_file))
        creation_time = os.stat(self.input_file).st_ctime
        nxs = self.nxs = Nexus(self.output_file, mode="w", creator="dahu")

        entry_grp = nxs.new_entry("entry", self.input.get("plugin_name", "dahu"),
                              title='BioSaxs multiframe integration',
                              force_time=get_isotime(creation_time))
        nxs.h5.attrs["default"] = entry_grp.name

        # Configuration
        cfg_grp = nxs.new_class(entry_grp, "configuration", "NXnote")
        cfg_grp.create_dataset("data", data=json.dumps(self.input, indent=2, separators=(",\r\n", ":\t")))
        cfg_grp.create_dataset("format", data="text/json")

        # Process 0: Measurement group
        measurement_grp = nxs.new_class(entry_grp, "0_measurement", "NXdata")
        measurement_grp.attrs["SILX_style"] = SAXS_STYLE
        # Instrument
        instrument_grp = nxs.new_instrument(entry_grp, "BM29")
        instrument_grp["name"] = "BioSaxs"
        source_grp = nxs.new_class(instrument_grp, "ESRF", "NXsource")
        source_grp["type"] = "Synchrotron X-ray source"
        source_grp["name"] = "European Synchrotron Radiation Facility"
        source_grp["probe"] = "X-ray"
        current = numpy.ascontiguousarray(self.input.get("storage_ring_current", []), dtype=numpy.float32)
        current_ds = source_grp.create_dataset("current",
                                               data=current)
        current_ds.attrs["units"] = "mA"
        current_ds.attrs["interpretation"] = "spectrum"

        # Sample: outsourced !
        create_nexus_sample(self.nxs, entry_grp, self.sample)

        monochromator_grp = nxs.new_class(instrument_grp, "multilayer", "NXmonochromator")
        wl_ds = monochromator_grp.create_dataset("wavelength", data=numpy.float32(self.ai.wavelength * 1e10))
        wl_ds.attrs["units"] = "Å"
        wl_ds.attrs["resolution"] = 0.014
        nrj_ds = monochromator_grp.create_dataset("energy", data=self.energy)
        nrj_ds.attrs["units"] = "keV"
        nrj_ds.attrs["resolution"] = 0.014

        detector_grp = nxs.new_class(instrument_grp, self.ai.detector.name, "NXdetector")
        dist_ds = detector_grp.create_dataset("distance", data=self.ai.dist)
        dist_ds.attrs["units"] = "m"
        xpix_ds = detector_grp.create_dataset("x_pixel_size", data=self.ai.pixel2)
        xpix_ds.attrs["units"] = "m"
        ypix_ds = detector_grp.create_dataset("y_pixel_size", data=self.ai.pixel1)
        ypix_ds.attrs["units"] = "m"
        f2d = self.ai.getFit2D()
        xbc_ds = detector_grp.create_dataset("beam_center_x", data=f2d["centerX"])
        xbc_ds.attrs["units"] = "pixel"
        ybc_ds = detector_grp.create_dataset("beam_center_y", data=f2d["centerY"])
        ybc_ds.attrs["units"] = "pixel"
        mask = self.ai.detector.mask
        mask_ds = detector_grp.create_dataset("pixel_mask", data=mask, **cmp_int)
        mask_ds.attrs["interpretation"] = "image"
        mask_ds.attrs["long_name"] = "Mask for invalid/hidden pixels"
        mask_ds.attrs["filename"] = self.input.get("mask_file")
        ct_ds = detector_grp.create_dataset("count_time", data=self.input.get("exposure_time"))
        ct_ds.attrs["units"] = "s"
        timestamps = self.input.get("timestamps")
        if not timestamps:
            nframes = self.input_frames.shape[0]
            start = time.mktime(time.strptime(self._start_time, "%Y-%m-%dT%H:%M:%SZ"))
            stop = time.mktime(time.strptime(self._end_time, "%Y-%m-%dT%H:%M:%SZ"))
            timestamps = numpy.linspace(start, stop, nframes)
        time_ds = detector_grp.create_dataset("timestamps",
                                              data=numpy.ascontiguousarray(timestamps, dtype=numpy.float64))
        time_ds.attrs["units"] = "s"
        time_ds.attrs["interpretation"] = "spectrum"
        frame_ds = detector_grp.create_dataset("frame_ids",
                                              data=numpy.ascontiguousarray(self.input.get("frame_ids", []), dtype=numpy.uint32))
        frame_ds.attrs["long_name"] = "Frame number"
        frame_ds.attrs["interpretation"] = "spectrum"
        if self.COPY_IMAGES:
            data = self.input_frames
            frames_ds = detector_grp.create_dataset("frames",
                                                     data=data,
                                                     chunks=(1,) + data.shape[-2:],
                                                     **cmp_int)
            frames_ds.attrs["interpretation"] = "image"
            measurement_grp["images"] = frames_ds
        else:  # use external links
            with Nexus(self.input_file, "r", timeout=self.timeout) as nxsr:
                entry = nxsr.get_entries()[0]
                if "measurement" in entry:
                    measurement = entry["measurement"]
                else:
                    self.log_error("No measurement in entry: %s of data_file: %s" % (entry, self.input_file))
                h5path = measurement["data"].name
            rel_path = os.path.relpath(os.path.abspath(self.input_file), os.path.dirname(os.path.abspath(self.output_file)))
            measurement_grp["images"] = detector_grp["frames"] = h5py.ExternalLink(rel_path, h5path)

        measurement_grp.attrs["signal"] = "images"

        diode_grp = nxs.new_class(instrument_grp, "beamstop_diode", "NXdetector")
        diode_ds = diode_grp.create_dataset("diode",
                                            data=numpy.ascontiguousarray(self.monitor_values, numpy.float32))
        diode_ds.attrs["interpretation"] = "spectrum"
        diode_ds.attrs["comment"] = "I1 = raw flux (I0) multiplied with the absorption of the sample"
        nf_ds = diode_grp.create_dataset("normalization_factor", data=self.normalization_factor)
        nf_ds.attrs["comment"] = "used to convert in abolute scattering"

        # few hard links
        measurement_grp["diode"] = diode_ds
        measurement_grp["timestamps"] = time_ds
        measurement_grp["ring_curent"] = current_ds

    # Process 1: pyFAI
        integration_grp = nxs.new_class(entry_grp, "1_integration", "NXprocess")
        integration_grp["sequence_index"] = 1
        integration_grp["program"] = "pyFAI"
        integration_grp["version"] = pyFAI.version
        integration_grp["date"] = get_isotime()
        cfg_grp = nxs.new_class(integration_grp, "configuration", "NXnote")
        cfg_grp.create_dataset("data", data=json.dumps(self.ai.get_config(), indent=2, separators=(",\r\n", ": ")))
        cfg_grp.create_dataset("format", data="text/json")
        cfg_grp.create_dataset("file_name", data=self.poni)
        pol_ds = cfg_grp.create_dataset("polarization_factor", data=polarization_factor)
        pol_ds.attrs["comment"] = "Between -1 and +1, 0 for circular"
        cfg_grp.create_dataset("integration_method", data=json.dumps(method.method._asdict()))
        integration_data = nxs.new_class(integration_grp, "results", "NXdata")
        integration_grp.attrs["title"] = str(self.sample)

    # Stage 1 processing: Integration frame per frame
        integrate1_results = self.process1_integration(self.input_frames)
        radial_unit, unit_name = str(self.unit).split("_", 1)
        q = numpy.ascontiguousarray(integrate1_results.radial, numpy.float32)
        I = numpy.ascontiguousarray(integrate1_results.intensity, dtype=numpy.float32)
        sigma = numpy.ascontiguousarray(integrate1_results.sigma, dtype=numpy.float32)

        self.to_memcached[radial_unit] = q
        self.to_memcached["I"] = I
        self.to_memcached["sigma"] = sigma

        q_ds = integration_data.create_dataset(radial_unit, data=q)
        q_ds.attrs["units"] = unit_name
        q_ds.attrs["long_name"] = "Scattering vector q (nm⁻¹)"

        int_ds = integration_data.create_dataset("I", data=I)
        std_ds = integration_data.create_dataset("errors", data=sigma)
        integration_data.attrs["signal"] = "I"
        integration_data.attrs["axes"] = [".", radial_unit]
        integration_data.attrs["SILX_style"] = SAXS_STYLE

        int_ds.attrs["interpretation"] = "spectrum"
        int_ds.attrs["units"] = "arbitrary"
        int_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"
        # int_ds.attrs["uncertainties"] = "errors" This does not work
        int_ds.attrs["scale"] = "log"
        std_ds.attrs["interpretation"] = "spectrum"

        hplc_data = nxs.new_class(integration_grp, "hplc", "NXdata")
        hplc_data.attrs["title"] = "Chromatogram"
        sum_ds = hplc_data.create_dataset("sum", data=numpy.ascontiguousarray(integrate1_results.intensity.sum(axis=-1), dtype=numpy.float32))
        sum_ds.attrs["interpretation"] = "spectrum"
        sum_ds.attrs["long_name"] = "Summed Intensity"
        hplc_data["frame_ids"] = frame_ds
        hplc_data.attrs["signal"] = "sum"
        hplc_data.attrs["axes"] = "frame_ids"

        if self.input.get("hplc_mode"):
            entry_grp.attrs["default"] = entry_grp.attrs["default"] = integration_grp.attrs["default"] = hplc_data.name
            self.log_warning("HPLC mode detected, stopping after frame per frame integration")
            return

        integration_grp.attrs["default"] = integration_data.name

    # Process 2: Freesas cormap
        cormap_grp = nxs.new_class(entry_grp, "2_correlation_mapping", "NXprocess")
        cormap_grp["sequence_index"] = 2
        cormap_grp["program"] = "freesas.cormap"
        cormap_grp["version"] = freesas.version
        cormap_grp["date"] = get_isotime()
        cormap_data = nxs.new_class(cormap_grp, "results", "NXdata")
        cormap_data.attrs["SILX_style"] = NORMAL_STYLE
        cfg_grp = nxs.new_class(cormap_grp, "configuration", "NXcollection")

        fidelity_abs = self.input.get("fidelity_abs", 0)
        fidelity_rel = self.input.get("fidelity_rel", 0)
        cfg_grp["fidelity_abs"] = fidelity_abs
        cfg_grp["fidelity_rel"] = fidelity_rel

    # Stage 2 processing
        cormap_results = self.process2_cormap(integrate1_results.intensity, fidelity_abs, fidelity_rel)
        cormap_data.attrs["signal"] = "probability"
        cormap_ds = cormap_data.create_dataset("probability", data=cormap_results.probability)
        cormap_ds.attrs["interpretation"] = "image"
        cormap_ds.attrs["long_name"] = "Probability to be the same"

        count_ds = cormap_data.create_dataset("count", data=cormap_results.count)
        count_ds.attrs["interpretation"] = "image"
        count_ds.attrs["long_name"] = "Longest sequence where curves do not cross each other"

        to_merge_ds = cormap_data.create_dataset("to_merge", data=numpy.arange(*cormap_results.tomerge, dtype=numpy.uint16))
        to_merge_ds.attrs["long_name"] = "Index of equivalent frames"
        cormap_grp.attrs["default"] = cormap_data.name
        if self.ispyb.url:
            self.to_pyarch["merged"] = cormap_results.tomerge

    # Process 3: time average and standard deviation
        average_grp = nxs.new_class(entry_grp, "3_time_average", "NXprocess")
        average_grp["sequence_index"] = 3
        average_grp["program"] = fully_qualified_name(self.__class__)
        average_grp["version"] = __version__
        average_data = nxs.new_class(average_grp, "results", "NXdata")
        average_data.attrs["SILX_style"] = SAXS_STYLE
        average_data.attrs["signal"] = "intensity_normed"

    # Stage 3 processing
        res3 = self.process3_average(cormap_results.tomerge)

        Iavg = numpy.ascontiguousarray(res3.average, dtype=numpy.float32)
        sigma_avg = numpy.ascontiguousarray(res3.deviation, dtype=numpy.float32)
        norm = numpy.ascontiguousarray(res3.normalization, dtype=numpy.float32)

        int_avg_ds = average_data.create_dataset("intensity_normed",
                                                  data=Iavg,
                                                  **cmp_float)
        int_avg_ds.attrs["interpretation"] = "image"
        int_avg_ds.attrs["formula"] = "sum_i(signal_i))/sum_i(normalization_i)"
        int_std_ds = average_data.create_dataset("intensity_std",
                                                   data=sigma_avg,
                                                   **cmp_float)
        int_std_ds.attrs["interpretation"] = "image"
        int_std_ds.attrs["formula"] = "sqrt(sum_i(variance_i)/sum_i(normalization_i))"
        int_std_ds.attrs["method"] = "Propagated error from weighted mean assuming poissonian behavour of every data-point"
        
        int_nrm_ds = average_data.create_dataset("normalization", data=norm)
        int_nrm_ds.attrs["formula"] = "sum_i(normalization_i))"
        
        average_grp.attrs["default"] = average_data.name

    # Process 4: Azimuthal integration of the time average image
        ai2_grp = nxs.new_class(entry_grp, "4_azimuthal_integration", "NXprocess")
        ai2_grp["sequence_index"] = 4
        ai2_grp["program"] = "pyFAI"
        ai2_grp["version"] = pyFAI.version
        ai2_grp["date"] = get_isotime()
        ai2_data = nxs.new_class(ai2_grp, "results", "NXdata")
        ai2_data.attrs["signal"] = "I"
        ai2_data.attrs["axes"] = radial_unit
        ai2_data.attrs["SILX_style"] = SAXS_STYLE
        ai2_data.attrs["title"] = str(self.sample)

        ai2_grp["configuration"] = integration_grp["configuration"]
        # ai2_grp["polarization_factor"] = integration_grp["polarization_factor"]
        # ai2_grp["integration_method"] = integration_grp["integration_method"]
        ai2_grp.attrs["default"] = ai2_data.name

    # Stage 4 processing
        intensity_std = res3.deviation
        if numexpr is None:
            variance = intensity_std * intensity_std
        else:
            variance = numexpr.evaluate("intensity_std**2")
        res2 = self.ai._integrate1d_ng(res3.average, self.npt,
                                       variance=variance,
                                       polarization_factor=polarization_factor,
                                       unit=self.unit,
                                       safe=False,
                                       method=method)
        if self.ispyb.url:
            self.to_pyarch["avg"] = res2

        _, self.to_memcached["I_avg"], self.to_memcached["sigma_avg"] = res2
        ai2_q_ds = ai2_data.create_dataset(radial_unit,
                                           data=numpy.ascontiguousarray(res2.radial, dtype=numpy.float32))
        ai2_q_ds.attrs["units"] = unit_name
        ai2_q_ds.attrs["long_name"] = "Scattering vector q (nm⁻¹)"

        ai2_int_ds = ai2_data.create_dataset("I", data=numpy.ascontiguousarray(res2.intensity, dtype=numpy.float32))
        ai2_std_ds = ai2_data.create_dataset("errors",
                                             data=numpy.ascontiguousarray(res2.sigma, dtype=numpy.float32))

        ai2_int_ds.attrs["interpretation"] = "spectrum"
        ai2_int_ds.attrs["units"] = "arbitrary"
        ai2_int_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"
        # ai2_int_ds.attrs["uncertainties"] = "errors" #this does not work
        ai2_std_ds.attrs["interpretation"] = "spectrum"
        ai2_int_ds.attrs["units"] = "arbitrary"
        # Finally declare the default entry and default dataset ...
        entry_grp.attrs["default"] = ai2_data.name

        # Export this to the output JSON
        # self.output["q"] = res2.radial
        # self.output["I"] = res2.intensity
        # self.output["std"] = res2.sigma

    def process1_integration(self, data):
        "First step of the processing, integrate all frames, return a IntegrationResult namedtuple"
        logger.debug("in process1_integration")
        intensity = numpy.empty((self.nb_frames, self.npt), dtype=numpy.float32)
        sigma = numpy.empty((self.nb_frames, self.npt), dtype=numpy.float32)
        idx = 0
        for i1, frame in zip(self.monitor_values, data):
            res = self.ai._integrate1d_ng(frame, self.npt,
                                          normalization_factor=i1 * self.scale_factor,
                                          error_model="poisson",
                                          polarization_factor=polarization_factor,
                                          unit=self.unit,
                                          safe=False,
                                          method=method)
            intensity[idx] = res.intensity
            sigma[idx] = res.sigma
            if self.ispyb.url:
                self.to_pyarch[idx] = res
            idx += 1
        if idx == 0:
            self.log_error(f"No frame iterated over in process1_integration! len(frames): {len(data)} len(monitor): {len(self.monitor_values)}", do_raise=False)
            radial = numpy.zeros(self.npt, dtype=numpy.float32)
        else:
            radial = res.radial
        return IntegrationResult(radial, intensity, sigma)

    def process2_cormap(self, curves, fidelity_abs, fidelity_rel):
        "Take the integrated data as input, returns a CormapResult namedtuple"
        logger.debug("in process2_cormap")
        count = numpy.empty((self.nb_frames, self.nb_frames), dtype=numpy.uint16)
        proba = numpy.empty((self.nb_frames, self.nb_frames), dtype=numpy.float32)
        for i in range(self.nb_frames):
            proba[i, i] = 1.0
            count[i, i] = 0
            for j in range(i):
                res = freesas.cormap.gof(curves[i], curves[j])
                proba[i, j] = proba[j, i] = res.P
                count[i, j] = count[j, i] = res.c
        tomerge = get_equivalent_frames(proba, fidelity_abs, fidelity_rel)
        return CormapResult(proba, count, tomerge)

    def process3_average(self, tomerge):
        "Average out the valid frames and return an AverageResult namedtuple"
        logger.debug("in process3_average")
        valid_slice = slice(*tomerge)
        mask = self.ai.detector.mask
        sum_data = (self.input_frames[valid_slice]).sum(axis=0)
        sum_norm = self.scale_factor * sum(self.monitor_values[valid_slice])    
        if numexpr is not None:
            # Numexpr is many-times faster than numpy when it comes to element-wise operations
            intensity_avg = numexpr.evaluate("where(mask==0, sum_data/sum_norm, 0.0)")
            intensity_std = numexpr.evaluate("where(mask==0, sqrt(sum_data)/sum_norm, 0.0)") # Assuming Poisson, no uncertainties on the diode
        else:
            with numpy.errstate(divide='ignore'):
                intensity_avg = sum_data / sum_norm
                intensity_std = numpy.sqrt(sum_data)/sum_norm
            wmask = numpy.where(mask)
            intensity_avg[wmask] = 0.0
            intensity_std[wmask] = 0.0
        return AverageResult(intensity_avg, intensity_std, sum_norm)

    def send_to_ispyb(self):
        if self.input.get("hplc_mode") == 0:
            if self.ispyb.url and parse_url(self.ispyb.url).host:
                ispyb = IspybConnector(*self.ispyb)
                ispyb.send_averaged(self.to_pyarch)
            else:
                self.log_warning(f"Not sending to ISPyB: no valid URL {self.ispyb.url}")

    def send_to_icat(self):
        #Some more metadata for iCat, as strings: 
        to_icat = copy.copy(self.to_pyarch)
        to_icat["experiment_type"] = "hplc" if self.input.get("hplc_mode") else "sample-changer"  
        to_icat["sample"] = self.sample
        to_icat["SAXS_maskFile"] = self.mask
        to_icat["SAXS_waveLength"] = str(self.ai.wavelength)
        to_icat["SAXS_normalisation"] = str(self.normalization_factor)
        to_icat["SAXS_diode_currents"] = str(self.monitor_values)
        to_icat["SAXS_numberFrames"] = str(self.nb_frames)
        to_icat["SAXS_timePerFrame"] = self.input.get("exposure_time", "?")
        to_icat["SAXS_detector_distance"] = str(self.ai.dist)
        to_icat["SAXS_pixelSizeX"] = str(self.ai.detector.pixel2)
        to_icat["SAXS_pixelSizeY"] = str(self.ai.detector.pixel1)
        f2d = self.ai.getFit2D()
        to_icat["SAXS_beam_center_x"] = str(f2d["centerX"])
        to_icat["SAXS_beam_center_y"] = str(f2d["centerY"])
        
        metadata = {"scanType": "integration"}
        return send_icat(sample=self.sample.name,
                         raw=os.path.dirname(os.path.dirname(os.path.abspath(self.input_file))),
                         path=os.path.dirname(os.path.abspath(self.output_file)),
                         data=to_icat, 
                         dataset = "integrate",
                         gallery=self.ispyb.gallery or os.path.join(os.path.dirname(os.path.abspath(self.output_file)), "gallery"), 
                         metadata=metadata)

    def send_to_memcached(self):
        "Send the content of self.to_memcached to the storage"
        dico={}
        key_base = self.output_file
        for k in sorted(self.to_memcached.keys(), key=lambda i:self.to_memcached[i].nbytes):
            key = f"{key_base}_{k}"
            dico[key] = json.dumps(self.to_memcached[k], cls=NumpyEncoder)
        return to_memcached(dico) 

