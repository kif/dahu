#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* SubtractBuffer: Search for the equivalence of buffers, average them and subtract from sample signal.  
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20/02/2020"
__status__ = "development"
__version__ = "0.1.0"

import os
import posixpath
import json
from collections import namedtuple
from dahu.plugin import Plugin
from dahu.factory import register
from dahu.cache import DataCache
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
from .common import Sample, Ispyb, get_equivalent_frames, cmp
from .nexus import Nexus, get_isotime

NexusJuice = namedtuple("NexusJuice", "filename h5path npt unit q I poni mask energy polarization signal2d error2d buffer concentration")
def read_nexus(self, filename):
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
        energy = mono_grp["energy"]
        img_grp = nxsr.get_class(entry_grp["3_time_average"], class_type="NXdata")
        image2d = img_grp["intensity_normed"]
        error2d = img_grp["intensity_std"]
        
    return NexusJuice(filename, h5path, npt, unit, q, I, poni, mask, energy, polarization, image2d, error2d)


class SubtractBuffer(Plugin):
    """Search for the equivalence of buffers, average them and subtract from sample signal.
    
        Typical JSON file:
    {
      "buffers_files" = ["buffer_001.h5", "buffer_002.h5"],
      "sample_file" = "sample.h5",
      "output_file" = "subtracted.h5"
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
        self.arrays = {}
        self.unit = None
        self.poni = None
        self.mask = None
        self.energy = None
        self.polarization_factor = None
           
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
        self.arrays = {}

    def process(self):
        Plugin.process(self)
        logger.debug("SubtractBuffer.process")
        self.create_nexus()
        
    def create_nexus(self):
        nxs = self.nxs = Nexus(self.output_file, mode="w")
        entry_grp = nxs.new_entry("entry", self.input.get("plugin_name", "dahu"), 
                              title='BioSaxs buffer subtraction', 
                              force_time=get_isotime())
        nxs.h5.attrs["default"] = entry_grp.name
        #Process 0: Measurement group
        input_grp = nxs.new_class(entry_grp, "0_measurement", "NXcollection")
        input_grp["sequence_index"] = 0

        with Nexus(self.sample_file, "r") as nxsr:
            entry = nxsr.get_entries()[0]
            h5path = entry.name
            nxdata_grp = nxsr.h5[entry.attrs["default"]]
            signal = nxdata_grp.attrs["signal"]
            axis = nxdata_grp.attrs["axes"]
            self.arrays["I"] = nxdata_grp[signal][()]
            self.arrays["q"] = nxdata_grp[axis][()]
            self.unit = pyFAI.units.to_unit(axis+"_"+nxdata_grp[axis].attrs["units"])
            self.npt = len(self.arrays["q"])
        input_grp["sample"] = h5py.ExternalLink(self.sample_file, h5path)

        for idx, buffer_file in enumerate(self.buffer_files):
                with Nexus(buffer_file, "r") as nxsr:
                    entry = nxsr.get_entries()[0]
                    h5path = entry.name
                input_grp["buffer_%i"%idx] = h5py.ExternalLink(buffer_file, h5path)     
        
        #Process 1: CorMap
            # Process 2: Freesas cormap 
        cormap_grp = nxs.new_class(entry_grp, "1_cormap", "NXprocess")
        cormap_grp["sequence_index"] = 1
        cormap_grp["program"] = "freesas"
        cormap_grp["version"] = freesas.version
        cormap_grp["date"] = get_isotime()
        cormap_data = nxs.new_class(cormap_grp, "results", "NXdata")
        
    # Stage 2 processing   
        cormap_results = self.process2_cormap(integrate1_results.intensity)
        
        cormap_data.attrs["signal"] = "probability"
        cormap_ds = cormap_data.create_dataset("probability", data=cormap_results.probability)
        cormap_ds.attrs["interpretation"] = "image"
        cormap_ds.attrs["long_name"] = "Probability to be the same"
        
        count_ds = cormap_data.create_dataset("count", data=cormap_results.count)
        count_ds.attrs["interpretation"] = "image"
        count_ds.attrs["long_name"] = "Longest sequence where curves do not cross each other"

        cormap_grp["to_merge"] = numpy.arange(*cormap_results.tomerge, dtype=numpy.uint16)
        cormap_grp["fidelity_abs"] = self.input.get("fidelity_abs", 0)
        cormap_grp["fidelity_rel"] = self.input.get("fidelity_rel", 0)
        cormap_grp.attrs["default"] = cormap_data.name

