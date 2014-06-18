#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement, print_function
__doc__ = """
Plugins for ID02:

Distortion correction
"""
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20140617"
__status__ = "development"
version = "0.2"

import sys, os, time, posixpath
import numpy
import logging
logger = logging.getLogger("dahu.id02")

from dahu.plugin import Plugin
from dahu.factory import register
from dahu.utils import get_isotime

import h5py
# silence non serious error messages, which are printed
# because we use h5py in a new thread (not in the main one)
# this is a bug seems to be fixed on newer version !!
# see this h5py issue 206: https://code.google.com/p/h5py/issues/detail?id=206
h5py._errors.silence_errors()

import pyFAI.distortion
if "TANGO_HOST" not in os.environ:
    raise RuntimeError("No TANGO_HOST defined")
import PyTango


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
            self.log_error("Replace output file %s"%output_hdf5_fn, False)
            os.unlink(output_hdf5_fn)
        out_hdf5 = h5py.File(output_hdf5_fn)
        self.output_ds = out_hdf5.create_dataset(output_hdf5_ds, shape, "float32",
                                                 chunks=(1,) + shape[1:],
                                                 maxshape=(None,) + shape[1:])
        spline = self.input.get("distortion_spline")
        if not os.path.isfile(spline):
            self.log_error("No spline file %s" % spline)
        detector = pyFAI.detectors.Detector(splineFile=spline)
        method = self.input.get("method", "lut")
        device = self.input.get("device", None)
        workgroup = self.input.get("workgroup", None)
        self.distortion = pyFAI.distortion.Distortion(detector, shape[-2:], method=method, device=device, workgroup=workgroup)
        self.distortion.calc_init()
        

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

@register
class Metadata(Plugin):
    """
    Plugin in charge of retrieving all metadata for ID02 and storing them into a HDF5 file
    """
    def __init__(self):
        Plugin.__init__(self)
        self.cycle = None
        self.c216 = None
        self.hdf5 = None
        self.hdf5_filename = None
        self.entry = None
        self.instrument = None
        self.group = None
        self.tfg_grp = None
        
    def setup(self, kwargs=None):
        """
        Structure of the input data:
        {
        "hdf5_filename":"/tmp/metadata.h5"
        "entry": "entry",
        "instrument":"id02"
        "c216":"id02/c216/0",
        }
        """
        Plugin.setup(self, kwargs)
        self.c216 = self.input.get("c216","id02/c216/0")
        self.cycle = self.input.get("cycle",1)
        if "hdf5_filename" not in self.input:
            self.log_error("hdf5_filename not in input")
        self.hdf5_filename = self.input.get("hdf5_filename")
        self.entry = self.input.get("entry","entry")
        self.instrument = self.input.get("instrument","id02")

    def process(self):
        
        self.create_hdf5()
        self.log_error("before error C216=%s TANGO_HOST=%s"%(self.c216,os.environ.get("TANGO_HOST")),False)
        c216ds = PyTango.DeviceProxy(str(self.c216))
        tfg = c216ds.GetTfuFrameTimes()
        cycle = c216ds.GetTfuCycles()
        self.tfg_grp["data"] = tfg
        self.tfg_grp["cycle"] = cycle
        self.tfg_grp["data"].attrs["interpretation"] = "scalar"
        #handle the case of multiple cycles: n*more frames, exposure time divided by n
        tmp=numpy.outer(numpy.ones(self.cycle),tfg).reshape(-1,2)/self.cycle
        self.tfg_grp["live_time"] = tmp[:,1]
        self.tfg_grp["live_time"].attrs["interpretation"] = "scalar"
        self.tfg_grp["dead_time"] = tmp[:,1]
        self.tfg_grp["dead_time"].attrs["interpretation"] = "scalar"
        

        
    def create_hdf5(self):
        """
        Create a HDF5 file and datastructure
        """
        try:
            self.hdf5 = h5py.File(self.hdf5_filename, 'a')
        except IOError as error:
            os.unlink(self.hdf5_filename)
            self.log_error("Unable to open %s: %s. Removing file and starting from scratch"%(self.hdf5_filename, error),False)
            self.hdf5 = h5py.File(self.hdf5_filename)


        if not self.entry.endswith("_"):
            self.entry += "_"
        entries = len([i.startswith(self.entry) for i in self.hdf5])
        self.entry = posixpath.join("","%s%04d" % (self.entry, entries))
        self.instrument = posixpath.join(self.entry, self.instrument)
        self.group = self.hdf5.require_group(self.instrument)
        self.group.parent.attrs["NX_class"] = "NXentry"
        self.group.attrs["NX_class"] = "NXinstrument"
        self.tfg_grp = self.hdf5.require_group(posixpath.join(self.instrument, "TFG"))
        self.tfg_grp["device"] = numpy.string_(self.c216)
        
        self.group.parent["title"] = numpy.string_("id02.metadata")
        self.group.parent["program"] = numpy.string_("Dahu")
        self.group.parent["start_time"] = numpy.string_(get_isotime())

    def teardown(self):
        if self.group:
            self.group.parent["end_time"] = numpy.string_(get_isotime())
        if self.hdf5:
            self.hdf5.close()
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
