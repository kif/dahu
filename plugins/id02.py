#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
logger = logging.getLogger("dahu.id02")

from dahu.plugin import Plugin
from dahu.factory import register

import h5py
import pyFAI, pyFAI.distortion

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
    def __init__(self):
        Plugin.__init__(self)
        #self.input_hdf5 = None
        self.input_ds = None
        self.output_ds = None
        self.distortion = None

    def setup(self, kwargs):
        """
        Structure of the input: {
            "input_hdf5_filename" : path,
            input_hdf5_dataset : path,
            output_hdf5_filename : path,
            output_hdf5_dataset : path,
            #distortion_dx : path,
            #distortion_dy : path,
            distortion_spline : path,
            ...
            }
        nota: dx is the row index, fast index, inside a line, small stride
        nota: dy is the line index, slow index, inside a column, large stride

        """
        self.input = kwargs
        input_hdf5_fn = kwargs.get("input_hdf5_filename", "")
        input_hdf5_ds = kwargs.get("input_hdf5_dataset", "")
        if not  os.path.isfile(input_hdf5_fn):
            self.log_error("No input HDF5 file %s" % input_hdf5_fn)
        in_hdf5 = h5py.File(input_hdf5_fn)
        if not input_hdf5_ds in in_hdf5:
            self.log_error("HDF5 file %s has no dataset" % (input_hdf5_fn, input_hdf5_ds))
        self.input_ds = in_hdf5[input_hdf5_ds]
        shape = self.input_ds.shape
        output_hdf5_fn = kwargs.get("output_hdf5_filename", "")
        output_hdf5_ds = kwargs.get("output_hdf5_dataset", "")
        if os.path.exists(output_hdf5_fn):
            os.unlink(output_hdf5_fn)
        out_hdf5 = h5py.File(output_hdf5_fn)
        self.output_ds = out_hdf5.create_dataset(output_hdf5_ds, shape, "float32",
                                                 chunks=(1,) + shape[1:],
                                                 maxshape=(None,) + shape[1:])
        spline = kwargs.get("distortion_spline")
        if not os.path.isfile(spline):
            self.log_error("No spline file %s" % spline)
        detector = pyFAI.detectors.Detector(splineFile=spline)
        self.distortion = pyFAI.distortion.Distortion(detector, shape[-2:])

    def process(self):
        """
        process every single frame in the HDF5 dataset
        """
        for i in range(self.input_ds.shape[0]):
            self.output_ds[i] = self.distortion.correct(self.input_ds[i])

    def teardown(self):
        if self.input_ds:
            self.input_ds.file.close()
        if self.output_ds:
            self.output_ds.file.close()
        self.distortion = None

if __name__ == "__main__":
    p = Distortion()
    t0 = time.time()
    p.setup({"input_hdf5_filename":"/nobackup/localtest/edf111.h5",
             "input_hdf5_dataset" : "/entry_0001/LImA_DATA/data",
             "output_hdf5_filename" : "/nobackup/localtest/test/edf111_proc.h5",
             "output_hdf5_dataset" : "/entry_0001/distortion_corrected/data",
             "distortion_spline": "/nobackup/localtest/distorsion.spline"
             })
    p.process()
    print("Processed %s frames in %.3fs" % (p.input_ds.shape[0], time.time() - t0))
    p.teardown()
