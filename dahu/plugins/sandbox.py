#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

from __future__ import with_statement, print_function

import PyTango
import h5py
import logging
import numpy
import os
import posixpath
import pyFAI
import pyFAI.distortion
import pyFAI.worker
import shutil
import sys
import time
import threading

from dahu.factory import register
from dahu.plugin import Plugin, plugin_from_function
from dahu.utils import get_isotime


__doc__ = """
Sandbox for Plugins for ID02:

* Distortion correction
"""
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "10/10/2014"
__status__ = "development"
version = "0.3"

logger = logging.getLogger("dahu.id02")


@register
class Distortion(Plugin):
    """This class performs the distortion correction

    later will be added:
    dark subtraction
    normlization for flatfield
    normalization for solid angle correction
    normalization for exposure time
    normalization for beam intensity

    Structure of the input: {
            "input_hdf5_filename" : path,
            "input_hdf5_dataset : path,
            "output_hdf5_filename" : path,
            "output_hdf5_dataset" : path,
            #distortion_dx : path,
            #distortion_dy : path,
            "distortion_spline" : path,
            ...
            method: "lut" or "csr",
            device: "cpu" or "gpu",
            }
        nota: dx is the row index, fast index, inside a line, small stride
        nota: dy is the line index, slow index, inside a column, large stride

    """
    def __init__(self):
        Plugin.__init__(self)
        self.input_ds = None
        self.output_ds = None
        self.distortion = None

    def setup(self, kwargs=None):
        """
        Structure of the input: {
            "input_hdf5_filename" : path,
            "input_hdf5_dataset : path,
            "output_hdf5_filename" : path,
            "output_hdf5_dataset" : path,
            #distortion_dx : path,
            #distortion_dy : path,
            "distortion_spline" : path,
            ...
            method: "lut" or "csr"
            device: "cpu" or "gpu"
            }
        nota: dx is the row index, fast index, inside a line, small stride
        nota: dy is the line index, slow index, inside a column, large stride
        """
        if kwargs is not None:
            self.input.update(kwargs)
        input_hdf5_fn = self.input.get("input_hdf5_filename", "")
        input_hdf5_ds = self.input.get("input_hdf5_dataset", "")
        if not  os.path.isfile(input_hdf5_fn):
            self.log_error("No input HDF5 file %s" % input_hdf5_fn)
        in_hdf5 = h5py.File(input_hdf5_fn)
        if not input_hdf5_ds in in_hdf5:
            self.log_error("HDF5 file %s has no dataset" % (input_hdf5_fn, input_hdf5_ds))
        self.input_ds = in_hdf5[input_hdf5_ds]
        shape = self.input_ds.shape
        output_hdf5_fn = self.input.get("output_hdf5_filename", "")
        output_hdf5_ds = self.input.get("output_hdf5_dataset", "")
        if os.path.exists(output_hdf5_fn):
            self.log_error("Replace output file %s" % output_hdf5_fn, False)
            os.unlink(output_hdf5_fn)
        out_hdf5 = h5py.File(output_hdf5_fn)
        self.output_ds = out_hdf5.create_dataset(output_hdf5_ds, shape, "float32",
                                                 chunks=(1,) + shape[1:],
                                                 maxshape=(None,) + shape[1:])
        spline = self.input.get("distortion_spline")
        if spline and os.path.isfile(spline):
            detector = pyFAI.detectors.Detector(splineFile=spline)
        else:
            self.log_error("No spline file %s" % spline, do_raise=False)
            detector = None
        method = self.input.get("method", "lut")
        device = self.input.get("device", None)
        workgroup = self.input.get("workgroup", 8)
        if detector:
            self.distortion = pyFAI.distortion.Distortion(detector, shape[-2:], method=method, device=device, workgroup=workgroup)
            self.distortion.calc_init()

    def process(self):
        """
        process every single frame in the HDF5 dataset
        """
        if self.distortion:
            for i in range(self.input_ds.shape[0]):
                self.output_ds[i] = self.distortion.correct(self.input_ds[i])
        else:
            for i in range(self.input_ds.shape[0]):
                self.output_ds[i] = self.input_ds[i]

    def teardown(self):
        if self.input_ds:
            self.input_ds.file.close()
        if self.output_ds:
            self.output_ds.file.close()
        self.distortion = None
        Plugin.teardown(self)


################################################################################
# Image filtering plugin
################################################################################
@register
class Filter(Plugin):
    """This plugin does filtering of a set of images taken in HDF5 and outputs a single image file (in TIF or EDF)

    @param image_file: HDF5 input file (mandatory)
    @param entry: entry in HDF5 input file (optional, default: last entry)
    @param output_file: name of the output file (mandatory)
    @param output_format: any format FabIO can write (default: edf)
    @param filter: "max", "min", "mean" or "median" (default: "mean")
    
    @param threshold: what is the upper limit? all pixel > max*(1-threshold) are discareded.
    @param minimum: minimum valid value or True
    @param maximum: maximum valid value
    @param darks: list of dark current images for subtraction
    @param flats: list of flat field images for division
    @param filter_: can be maximum, mean or median (default=mean)
    @param correct_flat_from_dark: shall the flat be re-corrected ?
    @param cutoff: keep all data where (I-center)/std < cutoff
    
    input = { "image_file": "/nobackup/lid02gpu11/FRELON/test_laurent_saxs_0000.h5",
              #"entry": "entry_0000"
              "output_format": "edf",
              "output_file": "/nobackup/lid02gpu12/dark.edf",
              "filter": "median",
              #"cutoff": 0
              #threshold:0,
              #"format": "edf",
              #"dark_file": filename,
              #"flat_file": filename,
#              "do_dark":false,
#              "do_flat":false,
              }  
    """
    
    def __init__(self):
        Plugin.__init__(self)
        self.images = None
        self.output_format = "edf"
        self.output_file = "toto.edf"
        self.filter = "mean" #Todo: average_percentil_20-80
        self.cutoff = None
        
    def setup(self):
        """
        Read unput parameters
        """
        #self.images = 
        if "output_format" in self.input:
            self.output_format = str(self.input["output_format"]).strip().lower()
        if "filter" in self.input:
            self.filter = str(self.input["filter"]).strip().lower()
        if "output_file" in self.input:
            self.output_file = self.input["output_file"]
           
    def process(self):
        flats = []
        darks = []    
        data = self.read_stack()
        if data is None:
            self.log_error("No dataset to average")
        pyFAI.utils.averageImages(data, filter_=self.filter, cutoff=self.cutoff,
                                threshold=0, format=self.output_format, output=self.output_file,
                                flats=flats, darks=darks)

    def read_stack(self):
        """
        read input dataset and populates self.images
        @return numpy array with the stack
        """
        nxs = pyFAI.io.Nexus(self.input["image_file"], "r")
        for entry in nxs.get_entries():
            for instrument in nxs.get_class(entry, class_type="NXinstrument"):
                for detector in nxs.get_class(instrument, class_type="NXdetector"):
                    for ds in nxs.get_data(detector):
                        print(ds)
                        return numpy.array(ds)

@register
class Peter(Plugin):
    """
    This plugin does all processing needed:
    - preprocess to generate the metadata entry
    - metadata retrieval using the C216 time frame generator
    - perform all transformation for all camera

    input = { TO be defined ....
            }
    """
    def __init__(self):
        Plugin.__init__(self)
        self.input_metadata = None
        self.plugins = {"metadata":Metadata()}

    def setup(self, kwargs=None):
        """
        see class documentation
        """
        Plugin.setup(self, kwargs)

    def process(self):
        """
        work a bit
        """
        self.input_metadata = preproc(**self.input)
        self.plugins["metadata"].input = self.input_metadata
        self.plugins["metadata"].setup()
        self.plugins["metadata"].process()
        self.plugins["metadata"].teardown()
        self.output_metadata =self.plugins["metadata"].o

    def teardown(self):
        Plugin.teardown(self)

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
