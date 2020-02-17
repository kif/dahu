#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* IntegrateMultiframe: perform the integration of many frames contained in a HDF5 file and average them  
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "17/02/2020"
__status__ = "development"
__version__ = "0.1.0"

import os
import json
from collections import namedtuple
from dahu.plugin import Plugin
from dahu.factory import register
from dahu.cache import DataCache
import logging
logger = logging.getLogger("bm29.integrate")
import numpy
try:
    import numexpr
except ImportError:
    logger.error("Numexpr is not installed, falling back on numpy's implementations")
    numexpr = None
import h5py
import fabio
import pyFAI, pyFAI.azimuthalIntegrator
from pyFAI.method_registry import IntegrationMethod
from hdf5plugin import Bitshuffle
import freesas, freesas.cormap
cmp = Bitshuffle()
#from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from .common import Sample, Ispyb, get_equivalent_frames
from .nexus import Nexus, get_isotime

KeyCache = namedtuple("KeyCache", "npt unit poni mask wavelength")
IntegrationResult = namedtuple("IntegrationResult", "radial intensity sigma")
CormapResult = namedtuple("CormapResult", "probability count tomerge")
AverageResult = namedtuple("AverageResult", "average deviation")


@register
class IntegrateMultiframe(Plugin):
    """perform the integration of many frames contained in a HDF5 file and average them
    
    Input parameters:
    :param poni_file: configuration of the geometry
    :param input_file: path for the HDF5 file
    :param 
    
    Typical JSON file:
    {
      "input_file": "/tmp/file1.h5",
      "output_file": "/tmp/file1.h5", # optional
      "max_frame": 1000,
      "frame_ids": [101, 102],
      "timestamps": [1580985678.47, 1580985678.58],
      "monitor_values": [1, 1.1],
      "storage_ring_current": [199.6, 199.5]
      "exposure_time": 0.1, 
      "normalisation_factor": 1.0,
      "poni_file": "/tmp/example.poni",
      "mask_file": "/tmp/mask.edf",
      "npt": 1000,
      "wavelength": 1e-10 #m
      "fidelity_abs": 0.1,
      "fidelity_rel": 0.5,
      "sample": {
        "name": "bsa",
        "description": "protein description like Bovine Serum Albumin",
        "buffer": "description of buffer, pH, ...",
        "concentration": 0,
        "hplc": "column name and chromatography conditions",
        "temperature": 20,
        "temperature_env": 20},  
      "ispyb": {
        "server": "http://ispyb.esrf.fr:1234",
        "login": "mx1234",
        "passwd": "secret",
        "pyarch_folder": "/data/pyarch/mx1234/1d", 
       } 
    }
    """
    cache = DataCache(5)
    COPY_IMAGES=False
    
    def __init__(self):
        Plugin.__init__(self)
        self.sample = None
        self.ispyb = None
        self.input_file = None
        self._input_frames = None
        self.output_file = None
        self.nxs = None 
        self.nb_frames = None
        self.ai = None
        self.npt = 1000
        self.unit = pyFAI.units.to_unit("q_nm^-1")
        self.unit_alt = pyFAI.units.to_unit("q_A^-1")
        self.polarization_factor = 0.9
        self.poni = self.mask = None
        self.wavelength = 1e-10 #A
        self.method = IntegrationMethod.select_method(1, "no", "csr", "opencl")[0]
        self.monitor_values = None
        self.normalization_factor = None

    def setup(self, kwargs=None):
        logger.debug("IntegrateMultiframe.setup")
        Plugin.setup(self, kwargs)
        
        self.ispyb = Ispyb._fromdict(self.input.get("ispyb", {}))
        self.sample = Sample._fromdict(self.input.get("sample", {}))
        
        self.input_file = self.input.get("input_file")
        if self.input_file is None:
            self.log_error("No input file provided", do_raise=True)
        self.output_file = self.input.get("output_file")
        if self.output_file is None:
            lst = list(os.path.splitext(self.input_file))
            lst.insert(1, "-integrate")
            self.output_file = "".join(lst)
            self.log_warning(f"No output file provided, using: {self.output_file}")
        self.nb_frames = len(self.input.get("frame_ids", []))
        self.npt = self.input.get("npt", self.npt)
        self.unit = self.input.get("unit", self.unit)
        self.poni = self.input.get("poni_file")
        if self.poni is None:
            self.log_error("No poni-file provided! aborting", do_raise=True)
        self.mask = self.input.get("mask_file")
        self.wavelength = self.input.get("wavelength", self.wavelength)
        self.monitor_values = numpy.array(self.input.get("monitor_values", 1), dtype=numpy.float64)
        self.normalization_factor = float(self.input.get("normalization_factor", 1))
        
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
        

    @property
    def input_frames(self):
        "For performance reasons, all frames are read in one bloc and cached, this returns a 3D numpy array"
        if self._input_frames is None:
            try:
                with Nexus(self.input_file, "r") as nxs:
                    entry = nxs.get_entries()[0]
                    measurement = nxs.get_class(entry, class_type="NXmeasurement")[0]
                    self._input_frames = measurement["data"][...]
            except Exception as err:
                self.log_error("%s: %s"%(type(err),str(err)), do_raise=True)
        return self._input_frames

    def process(self):
        "does the integration of a set of frames"
        logger.debug("IntegrateMultiframe.process")
        self.ai = self.create_integrator()
        self.create_nexus()
        #Send to ispyb

    def create_integrator(self):
        "use a cache if needed ... and return the integrator"
        key = KeyCache(self.npt, self.unit, self.poni, self.mask, self.wavelength)
        if key in self.cache:
            ai = self.cache[key]
        else:
            ai = pyFAI.load(self.poni)
            ai.wavelength = self.wavelength
            if self.mask:
                mask = numpy.logical_or(fabio.open(self.mask).data, ai.detector.mask).astype("int8")
                ai.detector.mask = mask
            self.cache[key] = ai
        return ai

    def create_nexus(self):
        creation_time = os.stat(self.input_file).st_ctime
        nxs = self.nxs = Nexus(self.output_file, mode="w")
        
        entry_grp = nxs.new_entry("entry", self.input.get("plugin_name", "dahu"), 
                              title='BioSaxs multiframe integration', 
                              force_time=get_isotime(creation_time))
        nxs.h5.attrs["default"] = entry_grp.name
        #Process 0: Measurement group
        measurement_grp = nxs.new_class(entry_grp, "0_measurement", "NXdata")
        
        #Instrument
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
    
        #Sample:
        sample_grp = nxs.new_class(entry_grp, self.sample.name, "NXsample")
        if self.sample.description is not None:
            sample_grp["description"] = self.sample.description
        if self.sample.concentration is not None:
            sample_grp["concentration"] = self.sample.concentration
            sample_grp["concentration"].attrs["units"] = "mg/mL"
        if self.sample.buffer is not None:
            sample_grp["buffer"] = self.sample.buffer
            sample_grp["buffer"].attrs["comment"] = "Buffer description"
        if self.sample.hplc:
            sample_grp["hplc"] =  self.sample.hplc
            sample_grp["hplc"].attrs["comment"] = "Conditions for HPLC experiment"
        if self.sample.temperature is not None:
            sample_grp["temperature"] = self.sample.temperature
            sample_grp["temperature"].attrs["units"] = "°C"
            sample_grp["temperature"].attrs["comment"] = "Exposure temperature"
        if self.sample.temperature_env is not None:
            sample_grp["temperature_env"] = self.sample.temperature_env
            sample_grp["temperature_env"].attrs["units"] = "°C"
            sample_grp["temperature_env"].attrs["comment"] = "Storage temperature"
        
        monochromator_grp = nxs.new_class(instrument_grp, "monochromator", "NXmonochromator")
        monochromator_grp["wavelength"] = self.ai.wavelength*1e10
        monochromator_grp["wavelength"].attrs["units"] = "Ångström" 
        monochromator_grp["wavelength"].attrs["resolution"] = 0.01
        monochromator_grp["wavelength"].attrs["monochromator"] = "multilayer"
        
        detector_grp = nxs.new_class(instrument_grp, str(self.ai.detector), "NXdetector")
        detector_grp["distance"] = self.ai.dist
        detector_grp["distance"].attrs["units"] = "m"
        detector_grp["x_pixel_size"] = self.ai.pixel2
        detector_grp["x_pixel_size"].attrs["units"] = "m"
        detector_grp["y_pixel_size"] = self.ai.pixel1
        detector_grp["y_pixel_size"].attrs["units"] ="m"
        f2d = self.ai.getFit2D()
        detector_grp["beam_center_x"] = f2d["centerX"]
        detector_grp["beam_center_x"].attrs["units"] = "pixel"
        detector_grp["beam_center_y"] = f2d["centerY"]
        detector_grp["beam_center_y"].attrs["units"] = "pixel"
        mask = self.ai.detector.mask
        mask_ds = detector_grp.create_dataset("pixel_mask",
                                               data=mask,
                                               **cmp)
        mask_ds.attrs["filename"] = self.input.get("mask_file")
        detector_grp["count_time"] = self.input.get("exposure_time")
        detector_grp["count_time"].attrs["units"] = "s"
        time_ds = detector_grp.create_dataset("timestamps",
                                              data=numpy.ascontiguousarray(self.input.get("timestamps", []), dtype=numpy.float64))
        time_ds.attrs["interpretation"] = "spectrum"
        if self.COPY_IMAGES:
            data = self.input_frames
            frames_ds = detector_grp.create_dataset("frames",
                                                     data=data,
                                                     chunks=(1,)+data.shape[-2:],
                                                     **cmp)
            frames_ds.attrs["interpretation"] = "image"
            measurement_grp["images"] = frames_ds
        else: #use external links
            with Nexus(self.input_file, "r") as nxs:
                entry = nxs.get_entries()[0]
                measurement = nxs.get_class(entry, class_type="NXmeasurement")[0]
                h5path = measurement["data"].name
            measurement_grp["images"] = detector_grp["frames"] = h5py.ExternalLink(self.input_file, h5path)            
            
        measurement_grp.attrs["signal"] = "images"
        
        diode_grp = nxs.new_class(instrument_grp, "beamstop_diode", "NXdetector")
        diode_ds = diode_grp.create_dataset("diode", 
                                            data = numpy.ascontiguousarray(self.monitor_values,numpy.float32))
        diode_ds.attrs["interpretation"] = "spectrum"
        diode_ds.attrs["comment"] = "I1 = raw flux (I0) multiplied with the absorption of the sample"
        diode_grp["normalization_factor"] = self.input.get("normalization_factor")
        diode_grp["normalization_factor"].attrs["comment"] = "used to convert in abolute scattering"
    

        
        measurement_grp["diode"] = diode_ds
 
        measurement_grp["timestamps"] = time_ds
        measurement_grp["ring_curent"] = current_ds
    
    # Process 1: pyFAI
        integration_grp = nxs.new_class(entry_grp, "1_integration", "NXprocess")
        integration_grp["sequence_index"] = 1
        integration_grp["program"] = "pyFAI"
        integration_grp["version"] = pyFAI.version
        integration_grp["date"] = get_isotime()
        integration_grp["configuration"] = json.dumps(self.ai.get_config())
        integration_grp["configuration"].attrs["format"] = "json"
        integration_grp["polarization_factor"] = self.polarization_factor
        integration_grp["polarization_factor"].attrs["comment"] = "Between -1 and +1, 0 for circular"
        integration_data = nxs.new_class(integration_grp, "results", "NXdata")
        integration_data.attrs["signal"] = "I"
        
    # Stage 1 processing: Integration frame per frame
        integrate1_results = self.process1_integration(self.input_frames)
        
        radial_unit, unit_name = str(self.unit).split("_", 1)
        q_ds = integration_data.create_dataset(radial_unit,
                                                data=numpy.ascontiguousarray(integrate1_results.radial, numpy.float32))
        q_ds.attrs["units"] = unit_name
#        radial_unit_alt, unit_name_alt = str(self.unit_alt).split("_", 1)
#        qalt_ds = integration_data.create_dataset(radial_unit_alt+,(self.npt,), 
#                                                  data=numpy.ascontiguousarray(integrate1_results.radial*self.unit_alt.scale, dtype=numpy.float32))
#        qalt_ds.attrs["units"] = unit_name_alt

        int_ds = integration_data.create_dataset("I", 
                                                 data=numpy.ascontiguousarray(integrate1_results.intensity, dtype=numpy.float32))
        std_ds = integration_data.create_dataset("errors", 
                                                 data=numpy.ascontiguousarray(integrate1_results.sigma, dtype=numpy.float32))
        integration_data.attrs["signal"] = "I"
        integration_data.attrs["axes"] = [".", radial_unit]
        int_ds.attrs["interpretation"] = "spectrum" 
        std_ds.attrs["interpretation"] = "spectrum"

        sum_ds = integration_data.create_dataset("sum", 
                                                 data=numpy.ascontiguousarray(integrate1_results.intensity.sum(axis=-1), dtype=numpy.float32))
        sum_ds.attrs["interpretation"] = "spectrum" 

        
    # Process 2: Freesas cormap 
        cormap_grp = nxs.new_class(entry_grp, "2_cormap", "NXprocess")
        cormap_grp["sequence_index"] = 2
        cormap_grp["program"] = "freesas"
        cormap_grp["version"] = freesas.version
        cormap_grp["date"] = get_isotime()
        cormap_data = nxs.new_class(cormap_grp, "results", "NXdata")
        
    # Stage 2 processing   
        cormap_results = self.process2_cormap(integrate1_results.intensity)
        
        cormap_data.attrs["signal"] = "probability"
        cormap_ds =  cormap_data.create_dataset("probability", 
                                                data=cormap_results.probability)
        cormap_ds.attrs["interpretation"] = "image"
        cormap_ds.attrs["units"] = "probability"
        
        count_ds = cormap_data.create_dataset("count", 
                                              data=cormap_results.count)
        count_ds.attrs["interpretation"] = "image"
        count_ds.attrs["units"] = "longest sequence where curves do not cross"

        cormap_grp["to_merge"] = numpy.arange(*cormap_results.tomerge, dtype=numpy.uint16)
        cormap_grp["fidelity_abs"] = self.input.get("fidelity_abs", 0)
        cormap_grp["fidelity_rel"] = self.input.get("fidelity_rel", 0)
        
    # Process 3: time average and standard deviation
        average_grp = nxs.new_class(entry_grp, "3_time_average", "NXprocess")
        average_grp["sequence_index"] = 3
        average_grp["program"] = "Weighted frame average"
        average_grp["version"] = __version__
        average_data = nxs.new_class(average_grp, "results", "NXdata")
        average_data.attrs["signal"] = "intensity_normed"
    
    # Stage 3 processing
        res3 = self.process3_average(cormap_results.tomerge)    
        int_avg_ds =  average_data.create_dataset("intensity_normed", 
                                                  data=numpy.ascontiguousarray(res3.average, dtype=numpy.float32),
                                                  **cmp)
        int_avg_ds.attrs["interpretation"] = "image"
        int_avg_ds.attrs["formula"] = "sum_i(signal_i))/sum_i(normalization_i)"
        int_std_ds =  average_data.create_dataset("intensity_std", 
                                                   data=numpy.ascontiguousarray(res3.deviation, dtype=numpy.float32),
                                                   **cmp)
        int_std_ds.attrs["interpretation"] = "image"    
        int_std_ds.attrs["formula"] = "sqrt(sum_i(variance_i))/sum(normalization_i)"
        int_std_ds.attrs["method"] = "Propagated error from weighted mean assuming poissonian behavour of every data-point"
    
    # Process 4: Azimuthal integration of the time average image
        ai2_grp = nxs.new_class(entry_grp, "4_azimuthal_integration", "NXprocess")
        ai2_grp["sequence_index"] = 4
        ai2_grp["program"] = "pyFAI"
        ai2_grp["version"] = pyFAI.version
        ai2_grp["date"] = get_isotime()
        ai2_data = nxs.new_class(ai2_grp, "results", "NXdata")
        ai2_data.attrs["signal"] = "I"
        ai2_grp["configuration"]=integration_grp["configuration"]
        ai2_grp["polarization_factor"] = integration_grp["polarization_factor"]

    # Stage 4 processing
        intensity_std = res3.deviation
        if numexpr is None:
            variance = intensity_std*intensity_std
        else:
            variance = numexpr.evaluate("intensity_std**2")
        res2 = self.ai._integrate1d_ng(res3.average, self.npt, 
                                       variance=variance,
                                       polarization_factor=self.polarization_factor,
                                       unit=self.unit,
                                       safe=False,
                                       method=self.method)

        ai2_q_ds = ai2_data.create_dataset(radial_unit,
                                           data=numpy.ascontiguousarray(res2.radial, dtype=numpy.float32))
        ai2_q_ds.attrs["units"] = unit_name
        ai2_int_ds = ai2_data.create_dataset("I", 
                                             data=numpy.ascontiguousarray(res2.intensity, dtype=numpy.float32))
        ai2_std_ds = ai2_data.create_dataset("errors", 
                                             data=numpy.ascontiguousarray(res2.sigma, dtype=numpy.float32))
        ai2_data.attrs["signal"] = "I"
        ai2_data.attrs["axes"] = [radial_unit]
        ai2_data.attrs["interpretation"] ="spectrum"     
    
        #Finally declare the default entry and default dataset ...
        entry_grp.attrs["default"] = ai2_data.name
    
    def process1_integration(self, data):
        "First step of the processing, integrate all frames, return a IntegrationResult namedtuple"
        logger.debug("in process1_integration")
        intensity = numpy.empty((self.nb_frames, self.npt), dtype=numpy.float32)
        sigma = numpy.empty((self.nb_frames, self.npt), dtype=numpy.float32)

        for idx, frame in enumerate(data):
            i1 = self.monitor_values[idx]/self.normalization_factor
            res = self.ai._integrate1d_ng(frame, self.npt, 
                                          normalization_factor=i1,
                                          error_model="poisson",
                                          polarization_factor=self.polarization_factor,
                                          unit=self.unit,
                                          safe=False,
                                          method=self.method)
            intensity[idx] = res.intensity
            sigma[idx] = res.sigma
        return IntegrationResult(res.radial, intensity, sigma)

    def process2_cormap(self, curves):
        "Take the integrated data as input, returns a CormapResult namedtuple"
        logger.debug("in process2_cormap")
        count = numpy.empty((self.nb_frames, self.nb_frames), dtype=numpy.uint16)
        proba = numpy.empty((self.nb_frames, self.nb_frames), dtype=numpy.float32)
        for i in range(self.nb_frames):
            proba[i, i] = 1.0
            count[i, i] = 0
            for j in range(i):
                res = freesas.cormap.gof(curves[i], curves[j])
                proba[i,j] = proba[j,i] = res.P
                count[i,j] = count[j,i] = res.c
        tomerge = get_equivalent_frames(proba, self.input.get("fidelity_abs", 0), self.input.get("fidelity_rel", 0))
        return CormapResult(proba, count, tomerge)
        
    def process3_average(self, tomerge):
        "Average out the valid frames and return an AverageResult namedtuple"
        logger.debug("in process3_average")
        valid_slice = slice(*tomerge)
        mask = self.ai.detector.mask
        sum_data = (self.input_frames[valid_slice]).sum(axis=0)
        sum_norm = (numpy.array(self.monitor_values)[valid_slice]).sum()/self.normalization_factor
        if numexpr is not None:
            #Numexpr is many-times faster than numpy when it comes to element-wise operations
            intensity_avg = numexpr.evaluate("where(mask==0, sum_data/sum_norm, 0.0)")
            intensity_std = numexpr.evaluate("where(mask==0, sqrt(sum_data)/sum_norm, 0.0)")
        else:            
            intensity_avg = sum_data/sum_norm
            intensity_std = numpy.sqrt(sum_data)/sum_norm
            wmask = numpy.where(mask)
            intensity_avg[wmask] = 0.0
            intensity_std[wmask] = 0.0
        return AverageResult(intensity_avg, intensity_std)
        