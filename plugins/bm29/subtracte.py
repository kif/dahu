#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* SubtractBuffer: Search for the equivalence of buffers, average them and subtract from sample signal.  
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "27/03/2020"
__status__ = "development"
__version__ = "0.1.0"

import os
import posixpath
import json
from collections import namedtuple
from dahu.plugin import Plugin
from dahu.factory import register
from dahu.cache import DataCache
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
import fabio
import pyFAI, pyFAI.azimuthalIntegrator
from pyFAI.method_registry import IntegrationMethod
import freesas, freesas.cormap
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from .common import Sample, Ispyb, get_equivalent_frames, cmp, get_integrator, KeyCache, polarization_factor, method
from .nexus import Nexus, get_isotime

NexusJuice = namedtuple("NexusJuice", "filename h5path npt unit q I poni mask energy polarization signal2d error2d buffer concentration")


class SubtractBuffer(Plugin):
    """Search for the equivalence of buffers, average them and subtract from sample signal.
    
        Typical JSON file:
    {
      "buffers_files" = ["buffer_001.h5", "buffer_002.h5"],
      "sample_file" = "sample.h5",
      "output_file" = "subtracted.h5"
      "fidelity" = 0.001,
      # TODO ... "wait_for" = ["job1", "job2"]
      
    } 
    """
    def __init__(self):
        Plugin.__init__(self)
        self.buffer_files = []
        self.sample_file = None
        self.nxs = None
        self.output_file = None
        self.ai = None
        self.npt = None
        self.poni = None
        self.mask = None
        self.energy = None
        self.sample_juice = None
        self.buffer_juices = []
           
    def setup(self, kwargs=None):
        logger.debug("SubtractBuffer.setup")
        Plugin.setup(self, kwargs)
        self.sample_file = self.input.get("sample_file")
        if self.sample_file is None:
            self.log_error("No sample file provided", do_raise=True)
        self.output_file = self.input.get("output_file")
        if self.output_file is None:
            lst = list(os.path.splitext(self.sample_file))
            lst.insert(1, "-sub")
            self.output_file = "".join(lst)
            self.log_warning(f"No output file provided, using: {self.output_file}")
        self.buffer_files = [os.path.abspath(fn) for fn in self.input.get("buffer_files", [])
                             if os.path.exists(fn)]
            
    def teardown(self):
        Plugin.teardown(self)
        logger.debug("SubtractBuffer.teardown")
        # export the output file location
        self.output["output_file"] = self.output_file
        if self.nxs is not None:
            self.nxs.close()
        if self.ai is not None:
            self.ai = None
        self.sample_juice = None
        self.buffer_juices = []

    def process(self):
        Plugin.process(self)
        logger.debug("SubtractBuffer.process")
        self.sample_juice = self.read_nexus(self.sample_file)
        self.create_nexus()
        
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
        if self.sample_juice.buffer != buffer_juice.buffer:
            self.log_warning(f"Sample {buffer_file} differs in buffer descsription, discarding")
            return
        if buffer_juice.concentration:
            self.log_warning(f"Sample {buffer_file} concentration not null, discarding")
            return
        return buffer_juice
    
    def create_nexus(self):
        nxs = self.nxs = Nexus(self.output_file, mode="w")
        entry_grp = nxs.new_entry("entry", self.input.get("plugin_name", "dahu"), 
                              title='BioSaxs buffer subtraction', 
                              force_time=get_isotime())
        nxs.h5.attrs["default"] = entry_grp.name
        #Process 0: Measurement group
        input_grp = nxs.new_class(entry_grp, "0_measurement", "NXcollection")
        input_grp["sequence_index"] = 0

        input_grp["sample"] = h5py.ExternalLink(self.sample_file, self.sample_juice.h5path)

        for idx, buffer_file in enumerate(self.buffer_files):
            buffer_juice = self.validate_buffer(buffer_file)
            if buffer_file is not None:
                input_grp["buffer_%i"%idx] = h5py.ExternalLink(buffer_file, buffer_juice.h5path)     
                self.buffer_juices.append(buffer_juice)
        
        #Process 1: CorMap
        cormap_grp = nxs.new_class(entry_grp, "1_cormap", "NXprocess")
        cormap_grp["sequence_index"] = 1
        cormap_grp["program"] = "freesas"
        cormap_grp["version"] = freesas.version
        cormap_grp["date"] = get_isotime()
        cormap_data = nxs.new_class(cormap_grp, "results", "NXdata")
        
        # cormap processing   
        nb_frames = len(self.buffer_juices)
        count = numpy.empty((nb_frames, nb_frames), dtype=numpy.uint16)
        proba = numpy.empty((nb_frames, nb_frames), dtype=numpy.float32)
        for i in range(nb_frames):
            proba[i, i] = 1.0
            count[i, i] = 0
            for j in range(i):
                res = freesas.cormap.gof(self.buffer_juices[i].I, self.buffer_juices[j].I)
                proba[i,j] = proba[j,i] = res.P
                count[i,j] = count[j,i] = res.c
        fidelity = self.input.get("fidelity", 0)
        tomerge = get_equivalent_frames(proba, fidelity, fidelity)

        cormap_data.attrs["signal"] = "probability"
        cormap_ds = cormap_data.create_dataset("probability", data=proba)
        cormap_ds.attrs["interpretation"] = "image"
        cormap_ds.attrs["long_name"] = "Probability to be the same"
        
        count_ds = cormap_data.create_dataset("count", data=count)
        count_ds.attrs["interpretation"] = "image"
        count_ds.attrs["long_name"] = "Longest sequence where curves do not cross each other"

        cormap_grp["to_merge"] = numpy.arange(*tomerge, dtype=numpy.uint16)
        cormap_grp["fidelity"] = fidelity
        cormap_grp.attrs["default"] = cormap_data.name
        
        #Process 2: Image processing: subtraction with standard deviation
        average_grp = nxs.new_class(entry_grp, "2_buffer_subtraction", "NXprocess")
        average_grp["sequence_index"] = 2
        average_grp["program"] = fully_qualified_name(self.__class__)
        average_grp["version"] = __version__
        average_data = nxs.new_class(average_grp, "results", "NXdata")
        average_data.attrs["signal"] = "intensity_normed"
        # Stage 2 

        #Nota: formula probably wrong ! Take into account the number of readings !
        #TODO implement those math using numexpr:
        #  avg = Σdata / Σnorm
        #  var = Σdata / (Σnorm)²
        #  Σnorm = avg / var
        #  Σdata = avg² / var
        # The propagate the error based on the number of frames in each buffer with quadratic summation
        # 
        buffer_average = numpy.mean([self.buffer_juices[i].signal2d for i in range(*tomerge)], axis=0)
        buffer_variance = numpy.sum([(self.buffer_juices[i].error2d)**2 for i in range(*tomerge)], axis=0) / (tomerge[1] - tomerge[0])**2
        sub_average = self.sample_juice.signal2d - buffer_average
        sub_variance = self.sample_juice.error2d**2 + buffer_variance
        sub_std = numpy.sqrt(sub_variance)
        
        
        int_avg_ds =  average_data.create_dataset("intensity_normed", 
                                                  data=numpy.ascontiguousarray(sub_average, dtype=numpy.float32),
                                                  **cmp)
        int_avg_ds.attrs["interpretation"] = "image"
        int_avg_ds.attrs["formula"] = "sample_signal - mean(buffer_signal_i)"
        int_std_ds =  average_data.create_dataset("intensity_std", 
                                                   data=numpy.ascontiguousarray(sub_std, dtype=numpy.float32),
                                                   **cmp)
        int_std_ds.attrs["interpretation"] = "image"    
        int_std_ds.attrs["formula"] = "sqrt( sample_variance + (sum(buffer_variance)/n_buffer**2 ))"
        int_std_ds.attrs["method"] = "quadratic sum of sample error and buffer errors"
        average_grp.attrs["default"] = average_data.name
        # Process 3: Azimuthal integration of the subtracted image
        ai2_grp = nxs.new_class(entry_grp, "3_azimuthal_integration", "NXprocess")
        ai2_grp["sequence_index"] = 3
        ai2_grp["program"] = "pyFAI"
        ai2_grp["version"] = pyFAI.version
        ai2_grp["date"] = get_isotime()
        key_cache = KeyCache(self.sample_juice.npt, self.sample_juice.unit, self.sample_juice.poni, self.sample_juice.mask, self.sample_juice.energy)
        ai = get_integrator(key_cache)
        ai2_grp.create_dataset("configuration", data=json.dumps(ai.get_config()))
        radial_unit, unit_name = str(key_cache.unit).split("_", 1)
        ai2_data = nxs.new_class(ai2_grp, "results", "NXdata")
        ai2_data.attrs["signal"] = "I"
        ai2_data.attrs["axes"] = radial_unit
        ai2_grp["polarization_factor"] = self.sample_juice.polarization
        ai2_grp.attrs["default"] = ai2_data.name

    # Stage 3 processing: azimuthal integration
        res2 = ai._integrate1d_ng(sub_average, key_cache.npt, 
                                  variance=sub_variance,
                                  polarization_factor=self.sample_juice.polarization,
                                  unit=key_cache.unit,
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
    
    #TODO:
    # Stage 4 processing: AutoRg --> available in FreeSAS
    # Stage 5 processing Kratky plot  --> copy from former EDNA
    # stage 6 Pair distribution function    --> Needs to implemented 
        
        #Finally declare the default entry and default dataset ...
        entry_grp.attrs["default"] = ai2_data.name


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
            npt = len(q)
            unit = pyFAI.units.to_unit(axis+"_"+nxdata_grp[axis].attrs["units"])
            integration_grp = nxdata_grp.parent
            poni = integration_grp["configuration"].attrs["poni_file"]
            polarization = integration_grp["polarization_factor"][()]
            instrument_grp = nxsr.get_class(entry_grp, class_type="NXinstrument")[0]
            detector_grp = nxsr.get_class(instrument_grp, class_type="NXdetector")[0]
            mask = detector_grp["pixel_mask"].attrs["filename"]
            mono_grp = nxsr.get_class(instrument_grp, class_type="NXmonochromator")[0]
            energy = mono_grp["energy"][()]
            img_grp = nxsr.get_class(entry_grp["3_time_average"], class_type="NXdata")[0]
            image2d = img_grp["intensity_normed"][()]
            error2d = img_grp["intensity_std"][()]
            sample_grp = nxsr.get_class(entry_grp, class_type="NXsample")[0]
            buffer = sample_grp["buffer"]
            concentration = sample_grp["concentration"]
        return NexusJuice(filename, h5path, npt, unit, q, I, poni, mask, energy, polarization, image2d, error2d, buffer, concentration)
