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
        self.mcs_grp = None
        
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
        self.read_c216()

    def create_hdf5(self):
        """
        Create a HDF5 file and datastructure
        """
        try:
            self.hdf5 = h5py.File(self.hdf5_filename, 'a')
        except IOError as error:
            os.unlink(self.hdf5_filename)
            self.log_error("Unable to open %s: %s. Removing file and starting from scratch" % (self.hdf5_filename, error), False)
            self.hdf5 = h5py.File(self.hdf5_filename)


        if not self.entry.endswith("_"):
            self.entry += "_"
        entries = len([i.startswith(self.entry) for i in self.hdf5])
        self.entry = posixpath.join("", "%s%04d" % (self.entry, entries))
        self.instrument = posixpath.join(self.entry, self.instrument)
        self.group = self.hdf5.require_group(self.instrument)
        self.group.parent.attrs["NX_class"] = "NXentry"
        self.group.attrs["NX_class"] = "NXinstrument"
        # TimeFrameGenerator
        self.tfg_grp = self.hdf5.require_group(posixpath.join(self.instrument, "TFG"))
        self.tfg_grp.attrs["NX_class"] = "NXcollection"
        self.tfg_grp["device"] = numpy.string_(self.c216)

        #MultiCounterScaler
        self.mcs_grp = self.hdf5.require_group(posixpath.join(self.instrument, "MCS"))
        self.mcs_grp.attrs["NX_class"] = "NXcollection"
        self.mcs_grp["device"] = numpy.string_(self.c216)

        #Factor
        HS32F = self.input.get("HS32F")
        if HS32F is not None:
            self.mcs_grp["HS32F"] = HS32F
        #Zero
        HS32Z = self.input.get("HS32Z")
        if HS32Z is not None:
            self.mcs_grp["HS32Z"] = HS32Z
        #Name
        HS32N = self.input.get("HS32N")
        if HS32N is not None:
            self.mcs_grp["HS32N"] = HS32N
        #Mode
        HS32M = self.input.get("HS32M")
        if HS32M is not None:
            self.mcs_grp["HS32M"] = HS32M
        if HS32N and HS32Z and HS32F and HS32M:
            self.mcs_grp.require_group("interpreted")
        self.group.parent["title"] = numpy.string_("id02.metadata")
        self.group.parent["program"] = numpy.string_("Dahu")
        self.group.parent["start_time"] = numpy.string_(get_isotime())

        
    def read_c216(self):
        """
        Manage the metadata coming from C216 Time Frame Generator
        """        
        c216ds = PyTango.DeviceProxy(str(self.c216))
        if c216ds.CompStatus("Tango::RUNNING") == "Tango::ON":
            msg = "C216 is running while reading counters ... possible issue"
            self._logging.append(msg)
            logger.warning(msg)
        
        raw_status = c216ds.GetCompleteConfigAndStatus()
        status = {"TFU_MODE":          raw_status[0],
                  "TFU_FRAME_MODE":    raw_status[1],
                  "START_MODE":        raw_status[2],
                  "FRAMES"     :       raw_status[14] // 2,
                  "CYCLES"      :      raw_status[17],
                  "CHANNELS"    :      raw_status[19],
                  "ACQ_BANK"       :   raw_status[65],
                  "TFU_CONFIG":        raw_status[:28],
                  "DIGI_IN":           raw_status[28:44],
                  "ANA_IN":            raw_status[44:60]}

        cycles = status["CYCLES"]
        frames = status["FRAMES"]
        tfg = c216ds.GetTfuFrameTimes()
        self.tfg_grp["frame_time"] = tfg
        self.tfg_grp["frame_time"].attrs["interpretation"] = "scalar"
        self.tfg_grp["cycles"] = cycles
        self.tfg_grp["frames"] = frames
        #handle the case of multiple cycles: n*more frames, exposure time always the same
        tmp = numpy.outer(numpy.ones(cycles),tfg)
        live_time =  tmp[1::2]
        self.tfg_grp["live_time"] = live_time
        self.tfg_grp["live_time"].attrs["interpretation"] = "scalar"
        self.tfg_grp["dead_time"] = tmp[0::2]
        self.tfg_grp["dead_time"].attrs["interpretation"] = "scalar"
        self.tfg_grp["delta_time"] = tmp.cumsum()[0::2]
        #raw scalers:
        raw_scalers = c216ds.ReadScalersForNLiveFrames([0, frames - 1])
        raw_scalers.shape = frames, -1
        counters = raw_scalers.shape[1]
        self.mcs_grp["HS32C"] = raw_scalers

        if "interpreted" in self.mcs_grp:
            modes = numpy.zeros(counters, dtype=numpy.int32)
            raw_modes = numpy.array(self.mcs_grp["HS32M"])
            modes[:raw_modes.size] = raw_modes
            values = numpy.zeros((frames,counters), dtype=numpy.float32)
            exptime = numpy.outer(tfg[1::2], numpy.ones(counters))
            zero = numpy.outer(numpy.ones(frames), numpy.array(self.mcs_grp["HS32Z"]))
            factor = numpy.outer(numpy.ones(frames), numpy.array(self.mcs_grp["HS32F"]))
            values_int = (raw_scalers-zero*exptime)*factor
            values_avg = (raw_scalers/exptime-zero)*factor           
            mask =(numpy.outer(numpy.ones(frames),modes)).astype(bool)
            values[mask] = values_int[mask]
            nmask = numpy.logical_not(mask)
            values[nmask] = values_avg[nmask]
            self.mcs_grp["HS32V"] = values.astype(numpy.float32)
            self.mcs_grp["HS32V"].attrs["interpretation"] = "scalar"
            for i, name in enumerate(self.mcs_grp["HS32N"]):
                fullname =  "interpreted/%s"%name
                self.mcs_grp[fullname]=values[:,i]
                self.mcs_grp[fullname].attrs["interpretation"] = "scalar"

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
