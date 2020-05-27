#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* IntegrateMultiframe: perform the integration of many frames contained in a HDF5 file and average them  
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "27/05/2020"
__status__ = "development"
__version__ = "0.2.0"

import os
import json
from collections import namedtuple
from dahu.plugin import Plugin
from dahu.factory import register
from dahu.utils import fully_qualified_name
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
import freesas, freesas.cormap

#from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from .common import Sample, Ispyb, get_equivalent_frames, get_integrator, KeyCache,\
                    method, polarization_factor,Nexus, get_isotime, cmp_float, cmp_int


IntegrationResult = namedtuple("IntegrationResult", "radial intensity sigma")
CormapResult = namedtuple("CormapResult", "probability count tomerge")
AverageResult = namedtuple("AverageResult", "average deviation")


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
      "energy": 12.0, #keV
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
        # self.polarization_factor = 0.9 --> constant
        self.poni = self.mask = None
        self.energy = None 
        #self.method = IntegrationMethod.select_method(1, "no", "csr", "opencl")[0] -> constant
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
        self.unit = pyFAI.units.to_unit(self.input.get("unit", self.unit))
        self.poni = self.input.get("poni_file")
        if self.poni is None:
            self.log_error("No poni-file provided! aborting", do_raise=True)
        self.mask = self.input.get("mask_file")
        self.energy = self.input.get("energy")
        if self.energy is None:
            self.log_error("No energy provided! aborting", do_raise=True)
        else:
            self.energy = numpy.float32(self.energy) #It is important to fix the datatype of the energy
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
                    if "measurement" in entry:
                        measurement = entry["measurement"]
                    else:
                        self.log_error("No measurement in entry: %s of data_file: %s" % (entry, self.input_file))
                    self._input_frames = measurement["data"][...]
            except Exception as err:
                self.log_error("%s: %s"%(type(err),str(err)), do_raise=True)
        return self._input_frames

    def process(self):
        "does the integration of a set of frames"
        logger.debug("IntegrateMultiframe.process")
        self.ai = get_integrator(KeyCache(self.npt, self.unit, self.poni, self.mask, self.energy))
        self.create_nexus()
        #Send to ispyb

    def create_nexus(self):
        creation_time = os.stat(self.input_file).st_ctime
        nxs = self.nxs = Nexus(self.output_file, mode="w", creator="dahu")
        
        entry_grp = nxs.new_entry("entry", self.input.get("plugin_name", "dahu"), 
                              title='BioSaxs multiframe integration', 
                              force_time=get_isotime(creation_time))
        nxs.h5.attrs["default"] = entry_grp.name
        
        # Configuration
        cfg_grp = nxs.new_class(entry_grp, "configuration", "NXnote")
        cfg_grp.create_dataset("data", data=json.dumps(self.input, indent=2, separators=(",\r\n", ":\t")))
        cfg_grp.create_dataset("format", data = "text/json")

        
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
            concentration_ds = sample_grp.create_dataset("concentration", data=self.sample.concentration)
            concentration_ds.attrs["units"] = "mg/mL"
        if self.sample.buffer is not None:
            buffer_ds = sample_grp.create_dataset("buffer", data=self.sample.buffer)
            buffer_ds.attrs["comment"] = "Buffer description"
        if self.sample.hplc:
            hplc_ds = sample_grp.create_dataset("hplc", data=self.sample.hplc)
            hplc_ds.attrs["comment"] = "Conditions for HPLC experiment"
        if self.sample.temperature is not None:
            tempe_ds = sample_grp.create_dataset("temperature", data=self.sample.temperature)
            tempe_ds.attrs["units"] = "°C"
            tempe_ds.attrs["comment"] = "Exposure temperature"
        if self.sample.temperature_env is not None:
            tempv_ds = sample_grp.create_dataset("temperature_env", data=self.sample.temperature_env)
            tempv_ds.attrs["units"] = "°C"
            tempv_ds.attrs["comment"] = "Storage temperature"
        
        monochromator_grp = nxs.new_class(instrument_grp, "multilayer", "NXmonochromator")
        wl_ds = monochromator_grp.create_dataset("wavelength", data=numpy.float32(self.ai.wavelength*1e10))
        wl_ds.attrs["units"] = "Å" 
        wl_ds.attrs["resolution"] = 0.014
        nrj_ds = monochromator_grp.create_dataset("energy", data=self.energy)
        nrj_ds.attrs["units"] = "keV" 
        nrj_ds.attrs["resolution"] = 0.014
        
        
        detector_grp = nxs.new_class(instrument_grp, str(self.ai.detector), "NXdetector")
        dist_ds = detector_grp.create_dataset("distance", data=self.ai.dist)
        dist_ds.attrs["units"] = "m"
        xpix_ds = detector_grp.create_dataset("x_pixel_size", data=self.ai.pixel2)
        xpix_ds.attrs["units"] = "m"
        ypix_ds = detector_grp.create_dataset("y_pixel_size", data= self.ai.pixel1)
        ypix_ds.attrs["units"] ="m"
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
        ct_ds = detector_grp.create_dataset("count_time",data=self.input.get("exposure_time"))
        ct_ds.attrs["units"] = "s"
        time_ds = detector_grp.create_dataset("timestamps",
                                              data=numpy.ascontiguousarray(self.input.get("timestamps", []), dtype=numpy.float64))
        time_ds.attrs["interpretation"] = "spectrum"
        if self.COPY_IMAGES:
            data = self.input_frames
            frames_ds = detector_grp.create_dataset("frames",
                                                     data=data,
                                                     chunks=(1,)+data.shape[-2:],
                                                     **cmp_int)
            frames_ds.attrs["interpretation"] = "image"
            measurement_grp["images"] = frames_ds
        else: #use external links
            with Nexus(self.input_file, "r") as nxsr:
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
                                            data = numpy.ascontiguousarray(self.monitor_values,numpy.float32))
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
        cfg_grp.create_dataset("format", data = "text/json")
        cfg_grp.create_dataset("file_name", data = self.poni)
        pol_ds = cfg_grp.create_dataset("polarization_factor", data=polarization_factor)
        pol_ds.attrs["comment"] = "Between -1 and +1, 0 for circular"
        cfg_grp.create_dataset("integration_method", data=json.dumps(method.method._asdict()))
        integration_data = nxs.new_class(integration_grp, "results", "NXdata")
        integration_grp.attrs["default"] = integration_data.name
        
        
    # Stage 1 processing: Integration frame per frame
        integrate1_results = self.process1_integration(self.input_frames)
        
        radial_unit, unit_name = str(self.unit).split("_", 1)
        q_ds = integration_data.create_dataset(radial_unit, data=numpy.ascontiguousarray(integrate1_results.radial, numpy.float32))
        q_ds.attrs["units"] = unit_name
        q_ds.attrs["long_name"] = "Scattering vector q (nm⁻¹)"

        int_ds = integration_data.create_dataset("I", data=numpy.ascontiguousarray(integrate1_results.intensity, dtype=numpy.float32))
        std_ds = integration_data.create_dataset("errors", data=numpy.ascontiguousarray(integrate1_results.sigma, dtype=numpy.float32))
        integration_data.attrs["signal"] = "I"
        integration_data.attrs["axes"] = [".", radial_unit]
        #integration_data.attrs["axes"] = numpy.array([".", radial_unit], dtype=h5py.string_dtype(encoding='utf-8'))
        
        int_ds.attrs["interpretation"] = "spectrum" 
        int_ds.attrs["units"] = "arbitrary"
        int_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"
        #int_ds.attrs["uncertainties"] = "errors" This does not work
        int_ds.attrs["scale"] = "log"
        std_ds.attrs["interpretation"] = "spectrum"

        sum_ds = integration_data.create_dataset("sum", data=numpy.ascontiguousarray(integrate1_results.intensity.sum(axis=-1), dtype=numpy.float32))
        sum_ds.attrs["interpretation"] = "spectrum" 

        
    # Process 2: Freesas cormap 
        cormap_grp = nxs.new_class(entry_grp, "2_correlation_mapping", "NXprocess")
        cormap_grp["sequence_index"] = 2
        cormap_grp["program"] = "freesas.cormap"
        cormap_grp["version"] = freesas.version
        cormap_grp["date"] = get_isotime()
        cormap_data = nxs.new_class(cormap_grp, "results", "NXdata")
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

        cormap_grp["to_merge"] = numpy.arange(*cormap_results.tomerge, dtype=numpy.uint16)
        cormap_grp.attrs["default"] = cormap_data.name
        
    # Process 3: time average and standard deviation
        average_grp = nxs.new_class(entry_grp, "3_time_average", "NXprocess")
        average_grp["sequence_index"] = 3
        average_grp["program"] = fully_qualified_name(self.__class__)
        average_grp["version"] = __version__
        average_data = nxs.new_class(average_grp, "results", "NXdata")
        average_data.attrs["signal"] = "intensity_normed"
        
    # Stage 3 processing
        res3 = self.process3_average(cormap_results.tomerge)    
        int_avg_ds =  average_data.create_dataset("intensity_normed", 
                                                  data=numpy.ascontiguousarray(res3.average, dtype=numpy.float32),
                                                  **cmp_float)
        int_avg_ds.attrs["interpretation"] = "image"
        int_avg_ds.attrs["formula"] = "sum_i(signal_i))/sum_i(normalization_i)"
        int_std_ds =  average_data.create_dataset("intensity_std", 
                                                   data=numpy.ascontiguousarray(res3.deviation, dtype=numpy.float32),
                                                   **cmp_float)
        int_std_ds.attrs["interpretation"] = "image"    
        int_std_ds.attrs["formula"] = "sqrt(sum_i(variance_i))/sum(normalization_i)"
        int_std_ds.attrs["method"] = "Propagated error from weighted mean assuming poissonian behavour of every data-point"
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

        ai2_grp["configuration"]=integration_grp["configuration"]
        ai2_grp.attrs["default"] = ai2_data.name

    # Stage 4 processing
        intensity_std = res3.deviation
        if numexpr is None:
            variance = intensity_std*intensity_std
        else:
            variance = numexpr.evaluate("intensity_std**2")
        res2 = self.ai._integrate1d_ng(res3.average, self.npt, 
                                       variance=variance,
                                       polarization_factor=polarization_factor,
                                       unit=self.unit,
                                       safe=False,
                                       method=method)

        ai2_q_ds = ai2_data.create_dataset(radial_unit,
                                           data=numpy.ascontiguousarray(res2.radial, dtype=numpy.float32))
        ai2_q_ds.attrs["units"] = unit_name
        ai2_q_ds.attrs["long_name"] = "Scattering vector q (nm⁻¹)"
        
        ai2_int_ds = ai2_data.create_dataset("I", data=numpy.ascontiguousarray(res2.intensity, dtype=numpy.float32))
        ai2_std_ds = ai2_data.create_dataset("errors", 
                                             data=numpy.ascontiguousarray(res2.sigma, dtype=numpy.float32))
        
        
        ai2_int_ds.attrs["interpretation"] ="spectrum"     
        ai2_int_ds.attrs["units"] = "arbitrary"
        ai2_int_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"
        #ai2_int_ds.attrs["uncertainties"] = "errors" #this does not work
        ai2_std_ds.attrs["interpretation"] ="spectrum"
        ai2_int_ds.attrs["units"] = "arbitrary"
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
                                          polarization_factor=polarization_factor,
                                          unit=self.unit,
                                          safe=False,
                                          method=method)
            intensity[idx] = res.intensity
            sigma[idx] = res.sigma
        return IntegrationResult(res.radial, intensity, sigma)

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
                proba[i,j] = proba[j,i] = res.P
                count[i,j] = count[j,i] = res.c
        tomerge = get_equivalent_frames(proba, fidelity_abs, fidelity_rel)
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
        