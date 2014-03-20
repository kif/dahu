#!/usr/bin/env python
# -*- coding: utf8 -*-

from __future__ import with_statement, print_function
__doc__ = """
Plugins for ID02:

Distortion correction and 
"""
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20140320"
__status__ = "development"
version = "0.1"

import sys, os, time
import logging
logger = loggin.getLogger("dahu.id02")

from dahu.plugin import Plugin 
from dahu.factory import register

import h5py
import pyFAI

@register
class Distortion(Plugin):
    """
    This class performs the distortion correction 
    
    later will be added:
    dark subtraction
    normlization for flatfield 
    normalization for solid angle correction
    normalization for exposure time
    normalization for beam intensity   
    """
    
    def setup(self, kwargs):
        """
        Structure of the input: {
            input_hdf5_filename = path,
            input_hdf5_dataset = path,
            output_hdf5_filename = path,
            output_hdf5_dataset = path,
            distortion_dx = path,
            distortion_dy = path,
            distortion_spline = path,
            ...
            }
        nota: dx is the row index, fast index, inside a line, small stride
        nota: dy is the line index, slow index, inside a column, large stride
        
        """
        
    def process(self):
         
    
    