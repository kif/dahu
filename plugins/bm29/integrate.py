#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* IntegrateMultiframe: perform the integration of many frames contained in a HDF5 file and average them  
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "07/02/2020"
__status__ = "development"
version = "0.0.1"

import os
import numpy
from dahu.plugin import Plugin
from dahu.factory import register
import logging
logger = logging.getLogger("plugin.bm29.integrate")


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
      "max_frame": 1000,
      "frame_ids": [101, 102],
      "time_stamps": [1580985678.47, 1580985678.58],
      "monitor_values": [1, 1.1],
      "exposure_time": 0.1,
      "normalisation_factor": 1.0,
      "poni_file": "/tmp/example.poni",
      "mask_file": "/tmp/mask.edf",
      "npt": 1000,
      "unit": "q_nm^-1",
      "fidelity_abs": 0.1,
      "fidelity_rel": 0.5
      "sample": {
        "code": "bsa",
        "comments": "protein name",
        "buffer": "description of buffer, pH, ...",
        "concentration": 0,
        "hplc": "column name and chromatography conditions",
        "storage_temp": 20,
        "exposure_temp": 20}, 
    }
    """
    def __init__(self):
        Plugin.__init__(self)

    def setup(self, kwargs):
        logger.debug("Integrate.setup")
        Plugin.setup(self, kwargs)