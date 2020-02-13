#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* IntegrateMultiframe: perform the integration of many frames contained in a HDF5 file and average them  
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "13/02/2020"
__status__ = "development"
version = "0.0.1"

import os
import numpy
import json
from collections import namedtuple
from dahu.plugin import Plugin
from dahu.factory import register
from dahu.cache import DataCache
import logging
logger = logging.getLogger("plugin.bm29.integrate")
import h5py
import fabio
import pyFAI, pyFAI.azimuthalIntegrator
from pyFAI.method_registry import IntegrationMethod
from hdf5plugin import Bitshuffle
import freesas
cmp = Bitshuffle()
#from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from .common import Sample, Ispyb, get_equivalent_frames
from .nexus import Nexus, get_isotime

KeyCache = namedtuple("KeyCache", "npt unit poni mask wavelength")
Integration1Result = namedtuple("Integration1Result", "radial signal error")
CormapResult = namedtuple("CormapResult", "probability count tomerge")

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
      "exposure_time": 0.1, 
      "normalisation_factor": 1.0,
      "poni_file": "/tmp/example.poni",
      "mask_file": "/tmp/mask.edf",
      "npt": 1000,
      "unit": "q_nm^-1", # optional
      "wavelength": 1 #Angstrom
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
        self.unit = "q_nm^-1"
        self.poni = self.mask = None
        self.wavelength = 1
        self.method = IntegrationMethod.select_method(1, "no", "csr", "opencl")[0]
        self.monitor_values = None

    def setup(self, kwargs):
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
        self.poni = self.input.get("poni")
        if self.poni is None:
            self.log_error("No poni-file provided! aborting", do_raise=True)
        self.mask = self.input.get("mask")
        self.wavelength = self.input.get("wavelength", self.wavelength)
        self.monitor_values = numpy.array(self.input.get("monitor_values", 1), dtype=numpy.float64)
        
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

    def create_integrator(self):
        "use a cache if needed ... and return the integrator"
        key = KeyCache(self.npt self.unit self.poni self.mask, self.wavelength)
        if key in self.cache:
            ai = self.cache[key]
        else:
            ai = pyFAI.load(self.poni)
            ai.wavelength = self.wavelength*1e-10
            if self.mask:
                mask = numpy.logical_or(fabio.open(self.mask).data, ai.detector.mask).astype("int8")
                ai.detector.mask = mask
            self.cache[key] = ai
        return ai

    def create_nexus(self):
        creation_time = os.stat(self.input_file).st_ctime
        nxs = self.nxs = Nexus(self.output_file, mode="w")
        
        entry_grp = nxs.new_entry("entry", 'dahu', 
                              title='BioSaxs Multiframe integration', 
                              force_time=get_isotime(creation_time))
        nxs.h5.attrs["default"] = entry_grp.name
        #Instrument
        instrument_grp = nxs.new_instrument(entry_grp, "BM29")
        instrument_grp["name"] = "BioSaxs"
        source_grp = nxs.new_class(instrument_grp, "ESRF", "NXsource")
        source_grp["type"] = "Synchrotron X-ray source"
        source_grp["name"] = "European Synchrotron Radiation Facility"
        source_grp["probe"] = "X-ray"
        current_ds = source_grp.require_dataset("current", 
                                               (self.nb_frames,), 
                                                dtype=numpy.float32)
        current_ds.attrs["unit"] = "mA"
    
        #Sample:
        sample_grp = nxs.new_class(entry_grp, self.sample.name, "NXsample")
        sample_grp["description"] = self.sample.description
        sample_grp["concentration"] = self.sample.concentration
        sample_grp["concentration"].attrs["unit"] = "mg/mL"
        sample_grp["buffer"] = self.sample.buffer
        sample_grp["buffer"].attrs["comment"] = "Buffer description"
        sample_grp["hplc"] =  self.sample.hplc
        sample_grp["hplc"].attrs["comment"] = "Conditions for HPLC experiment"
        sample_grp["temperature_env"] = self.sample.temperature_env
        sample_grp["temperature_env"].attrs["unit"] = "°C"
        sample_grp["temperature_env"].attrs["comment"] = "Storage temperature"
        sample_grp["temperature"] = self.sample.temperature
        sample_grp["temperature"].attrs["unit"] = "°C"
        sample_grp["temperature"].attrs["comment"] = "Exposure temperature"
        
        monochromator_grp = nxs.new_class(instrument_grp, "monochromator", "NXmonochromator")
        monochromator_grp["wavelength"] = self.ai.wavelength*1e10
        monochromator_grp["wavelength"].attrs["unit"] = "Ångström" 
        monochromator_grp["wavelength"].attrs["resolution"] = 0.01
        monochromator_grp["wavelength"].attrs["monochromator"] = "multilayer"
        
        #create links
        #instrument_grp["wavelength"] = monochromator_grp["wavelength"]
        #instrument_grp["normalization_factor"] = normalization_factor
    
        detector_grp = nxs.new_class(instrument_grp, str(self.ai.detector), "NXdetector")
        detector_grp["distance"] = self.ai.dist
        detector_grp["distance"].attrs["unit"] = "m"
        detector_grp["x_pixel_size"] = self.ai.pixel2
        detector_grp["x_pixel_size"].attrs["unit"] = "m"
        detector_grp["y_pixel_size"] = self.ai.pixel1
        detector_grp["y_pixel_size"].attrs["unit"] ="m"
        f2d = self.ai.getFit2D()
        detector_grp["beam_center_x"] = f2d["centerX"]
        detector_grp["beam_center_x"].attrs["unit"] = "pixel"
        detector_grp["beam_center_y"] = f2d["centerY"]
        detector_grp["beam_center_y"].attrs["units"] = "pixel"
        mask = self.ai.detector.mask
        mask_ds = detector_grp.require_dataset("pixel_mask",
                                               mask.shape,
                                               dtype=mask.dtype,
                                               chunks=mask.shape,
                                               **cmp)
        mask_ds[...] = mask
        mask_ds.attrs["filename"] = self.input.get("mask_file")
        detector_grp["count_time"] = self.input.get("exposure_time")
        detector_grp["count_time"].attrs["unit"] = "s"
        time_ds = detector_grp.require_dataset("timestamps",
                                               data=self.input.get("timestamps") 
                                               (self.nb_frames,), 
                                               dtype=numpy.float64)
        if self.COPY_IMAGES:
            data = self.input_frames
            frames_ds = detector_grp.require_dataset("frames",
                                                     data.shape,
                                                     dtype=data.dtype,
                                                     chunks=(1,)+data.shape[-2:],
                                                     data=data
                                                     **cmp)
            frames_ds.attrs["interpretation"] = "image"
        else: #use external links
            with Nexus(self.input_file, "r") as nxs:
                entry = nxs.get_entries()[0]
                measurement = nxs.get_class(entry, class_type="NXmeasurement")[0]
                h5path = measurement["data"].name
            detector_grp["frames"] = h5py.ExternalLink(self.input_file, h5path)
            
        
        diode_grp = nxs.new_class(instrument_grp, "beamstop_diode", "NXdetector")
        diode_ds = diode_grp.create_dataset("diode", 
                                            (self.nb_frames,),
                                            data = self.monitor_values,
                                            dtype=numpy.float32)
        diode_grp["normalization_factor"] = self.input.get("normalization_factor")
    
        # Process 1: pyFAI
        integration_grp = nxs.new_class(entry_grp, "1_integration", "NXprocess")
        integration_grp["sequence_index"] = 1
        integration_grp["program"] = "pyFAI"
        integration_grp["version"] = pyFAI.version
        integration_grp["date"] = get_isotime()
        integration_data = nxs.new_class(integration_grp, "results", "NXdata")
        integration_data.attrs["signal"] = "I"
        
        # Stage 1 processing: Integration frame per frame
        #TODO
        
        radial_unit, unit_name = str(self.unit).split("_", 1)
        q_ds = integration_data.require_dataset(radial_unit,(self.npt,), dtype=numpy.float32)
        q_ds.attrs["unit"] = unit_name
        int_ds = integration_data.require_dataset("I", (self.nb_frames, self.npt), dtype=numpy.float32)
        std_ds = integration_data.require_dataset("errors", (self.nb_frames, self.npt), dtype=numpy.float32)
        integration_data.attrs["signal"] = "I"
        integration_data.attrs["axes"] = [".", radial_unit]
        int_ds.attrs["interpretation"] ="spectrum" 
        std_ds.attrs["interpretation"] ="spectrum" 
        
        
         
        # Process 2: Freesas: 
        cormap_grp = nxs.new_class(entry_grp, "2_cormap", "NXprocess")
        cormap_grp["sequence_index"] = 2
        cormap_grp["program"] = "freesas"
        cormap_grp["version"] = freesas.version
        cormap_grp["date"] = get_isotime()
        cormap_data = nxs.new_class(cormap_grp, "results", "NXdata")
        
        cormap_results = self.process2_cormap(curves)
        
        cormap_data.attrs["signal"] = "probability"
        cormap_ds =  cormap_data.require_dataset("probability", 
                                                 shape=(self.nb_frames, self.nb_frames), 
                                                 dtype=numpy.float32)
        cormap_ds.attrs["interpretation"] = "image"
        cormap_ds.attrs["unit"] = "probability"
        cormap_ds[...] = cormap_results.probability
        
        count_ds =  cormap_data.require_dataset("count", 
                                                shape=(self.nb_frames, self.nb_frames),
                                                dtype=numpy.uint8)
        count_ds.attrs["interpretation"] = "image"
        count_ds.attrs["unit"] = "number of data-point"
        count_ds[...] = cormap_results.count

        cormap_data["to_merge"] = numpy.arange(cormap_results.tomerge)
        
        # process 3: time average and standard deviation
        average_grp = nxs.new_class(entry_grp, "3_time_average", "NXprocess")
        average_grp["sequence_index"] = 3
        average_grp["program"] = "experimental numpy implementation"
        average_grp["version"] = "0.0"
        average_data = nxs.new_class(average_grp, "results", "NXdata")
        average_data.attrs["signal"] = "intensity_normed"
    
        int_avg_ds =  average_data.require_dataset("intensity_normed", 
                                                  shape=mask.shape, 
                                                  dtype=numpy.float32,
                                                  **cmp
                                                  )
        int_avg_ds.attrs["interpretation"] = "image"
    
        int_std_ds =  average_data.require_dataset("intensity_error", 
                                                  shape=mask.shape, 
                                                  dtype=numpy.float32,
                                                  **cmp
                                                  )
        int_std_ds.attrs["interpretation"] = "image"    
        int_std_ds.attrs["method"] = "standard deviation over time-serie"
    
        # process 4: Azimuthal integration of the time average image
        ai2_grp = nxs.new_class(entry_grp, "4_azimuthal_integration", "NXprocess")
        ai2_grp["sequence_index"] = 4
        ai2_grp["program"] = "experimental pyFAI implementation"
        ai2_grp["version"] = "0.16"
        ai2_data = nxs.new_class(ai2_grp, "results", "NXdata")
        ai2_data.attrs["signal"] = "I"
#        ai2_data.attrs["signal"] = "I"
        ai2_q_ds = ai2_data.require_dataset(radial_unit,(self.npt,), dtype=numpy.float32)
        ai2_q_ds.attrs["unit"] = "nm^-1"
        ai2_int_ds = ai2_data.require_dataset("I", (self.npt,), dtype=numpy.float32)
        ai2_std_ds = ai2_data.require_dataset("errors", (self.npt,), dtype=numpy.float32)
        ai2_data.attrs["signal"] = "I"
        ai2_data.attrs["axes"] = [radial_unit]
        ai2_data.attrs["interpretation"] ="spectrum" 
    
        #Measurement group ...
        measurement_grp = nxs.new_class(entry_grp, "0_measurement", "NXdata")
        measurement_grp.attrs["signal"] = "images"
        measurement_grp["diode"] = diode_ds
        measurement_grp["images"] = frames_ds
        
        measurement_grp["timestamps"] = time_ds
        measurement_grp["ring_curent"] = current_ds
    
        integration_grp["configuration"] = json.dumps(self.ai.get_config())
        integration_grp["configuration"].attrs["format"] = "json"
       
        # Stage 1: Integration 
        
        q_ds[:] = res[0]
        int_ds = res[1]
        std_ds = res[2]
        integrated = res[1]
        
        
        integrated = numpy.zeros((self.nb_frames, self.npt), dtype=numpy.float32)
        for frame in 
        for x in xml:
            e = XmlParser(x)
            fimg = fabio.open(x[:-4]+".edf")
            idx = int(fimg.header.get("acq_frame_nb", 0))
            time_ds[idx] = float(fimg.header.get("time_of_day", 0))
            current_ds[idx] = e.extract("machineCurrent", float)
            print("set frame %i"%idx)
            frames_ds[idx] = fimg.data
            i1 = e.extract("beamStopDiode", float)
            diode_ds[idx] = e.extract("beamStopDiode", float)
            print("integrate frame %i"%idx)
            res = ai.integrate1d(fimg.data, npt, 
                                 normalization_factor=i1/normalization_factor,
                                 error_model="poisson",
                                 polarization_factor=polarization_factor,
                                 unit=unit,
                                 safe=False,
                                 method=method,
                                 )
            int_ds[idx] = res[1]
            std_ds[idx] = res[2]
            integrated[idx] = res[1]
        q_ds[:] = res[0]
        print(etree.tounicode(e.e))
        json.dumps(ai.get_config())
        
        
        #Processing for cormap
        
        for i in range(nframes):
            cormap_ds[i,i] = 1.0
            count_ds[i,i] = 0.0
            for j in range(i):
                res = cormap.gof(int_ds[i], int_ds[j])
                cormap_ds[i,j] = cormap_ds[j,i] = res.P
                count_ds[i,j] = count_ds[j,i] = res.c
                
        #processing for time_averaging
    
        
        #let's assume all frames are valid:
        valid = [ True ]*nframes
        
        frames = numpy.zeros((nframes,)+mask.shape)
        frames_norm = numpy.zeros_like(frames)
        i1 = numpy.zeros(nframes)
        idx = 0
        for i in range(nframes):
            if valid[i]:
                bsd = diode_ds[i][...]
                if bsd!=0.0:
                    frames_norm[idx, ...] = frames_ds[idx] * (normalization_factor/ bsd)
                    idx+=1
                    
        #cut the trailing part of the array if there are missing values 
        frames_norm = frames_norm[:idx]
        int_avg_ds[...] = numpy.mean(frames_norm, axis = 0)
        int_std_ds[...] = numpy.std(frames_norm, axis = 0)
        
        # Final azimuthal integration:
        
        polarization = ai.polarization(mask.shape, polarization_factor, with_checksum=False)
        solid_angle_correction = ai.solidAngleArray(shape=mask.shape, absolute=False)
        res1 = ai.integrate1d(polarization*solid_angle_correction, npt,
                               unit=unit,
                               method=method,
                               correctSolidAngle=False)
        denom = res1[1]
        ai2_q_ds[...] = res1[0]
        res2 = ai.integrate1d(int_avg_ds[...], npt,
                               unit=unit,
                               method=method,
                               correctSolidAngle=False)
        res3 = ai.integrate1d(int_std_ds[...]**2, npt,
                               unit=unit,
                               method=method,
                               correctSolidAngle=False)
        ai2_int_ds[...] = res2.sum/res1.sum
        ai2_std_ds[...] = numpy.sqrt(res3.sum/res3.count)/res1.sum
        
        #Finally declare the default entry and default dataset ...
        entry_grp.attrs["default"] = ai2_data.name
    
    def process1_integration(self):
        "First step of the processing, integrate all frames, return the Integration1Result"
        logger.debug("in process1_integration")
        signal = numpy.zeros((self.nb_frames, self.npt), dtype=numpy.float32)
        errors = numpy.zeros((self.nb_frames, self.npt), dtype=numpy.float32)

        for idx, frame in enumerate(self.input_frames):
            i1 = self.monitor_values[idx]/self.normalization_factor
            res = self.ai._integrate1d_ng(frame, self.npt, 
                                          normalization_factor=i1,
                                          error_model="poisson",
                                          polarization_factor=self.polarization_factor,
                                          unit=self.unit,
                                          safe=False,
                                          method=self.method)
            signal[idx] = res[1]
            errors[idx] = res[2]
        return Integration1Result(res.radial, signal, errors)

    def process2_cormap(self, curves):
        "Take the integrated data as input, returns a namedtuple"
        logger.debug("in process2_cormap")
        count = numpy.empty((self.nb_frames, self.nb_frames), dtype=numpy.uint8)
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
        