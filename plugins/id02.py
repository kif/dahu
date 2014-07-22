#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
from dahu.plugin import Plugin
from dahu.utils import get_isotime


__doc__ = """
Plugins for ID02:

* Distortion correction
* Metadata saving (C216)
* single detector prcessing
"""
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "22/07/2014"
__status__ = "development"
version = "0.3"

logger = logging.getLogger("dahu.id02")


# silence non serious error messages, which are printed
# because we use h5py in a new thread (not in the main one)
# this is a bug seems to be fixed on newer version !!
# see this h5py issue 206: https://code.google.com/p/h5py/issues/detail?id=206
try:
    h5py._errors.silence_errors()
except:
    pass

if "TANGO_HOST" not in os.environ:
    raise RuntimeError("No TANGO_HOST defined")


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

@register
class Metadata(Plugin):
    """
    Plugin in charge of retrieving all metadata for ID02 and storing them into a HDF5 file
    
    TODO: rewrite using Nexus class from pyFAI: shorter code

    NOTA: pin number are 1-based (I0, I1, time)

    Structure of the input data:
input = {
        "hdf5_filename":"/nobackup/lid02gpu11/metadata/test.h5",
        "entry": "entry",
        "instrument":"id02",
        "c216":"id02/c216/0",
        "HS32F": [1e-06, 1, 7763480, 8176290, 342239, 967341, 5541980, 1739160, 2753.61, 1351920, 140000000, 16719500, 1, 0.000995868, 0.000995868, 1],
        "HS32Z": [0, 0, 383.55, 126.4, 6582.1, 6973.6, 162.95, 0, 221.2, 207.8, 315.26, 145.11, 323.76, 170, 170, 228.4],
        "HS32N": ["TIME", "AUX1", "PIN41", "PIN42", "PIN5", "PIN6", "PIN7", "PIN8", "PIN1", "PIN2", "PIN3", "PIN4", "AUX2", "THC1", "THC2", "PRESS"],
        "HSI0": 12,
        "HSI1": 7,
        "HSTime": 1,
        "HMStartEpoch": 1405087717.12159,
        "HMStartTime": "2014-07-11T16:08:37.121591+0200",
        "info": {"DetectorInfo":"VxH:detbox=14952.235x0.000x1.000,dettab=-62.000x-245.000",
                 "ExperimentInfo":"0",
                 "MachineInfo": "Ie=183.27mA,u35u=100.000mm/0.001mm,u21m=100.000mm/0.000mm,u21d=100.000mm/-0.000mm",
                 "MirrorInfo": "rz=-3.600mrad,ty=0.455mm,ty1=2.075mm,ty2=-1.165mm,tz1=-0.030mm,tz2=-0.090mm,mir2rz=2.000mrad",
                 "OpticsInfo": "egy=12460.0eV,theta=9.132deg,hgt=11.7mm,tilt=4.440deg,tra=1.963mm",
                 "ProposalInfo": 0,
                 "StationInfo": "ID02"
                 }
        }

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
        see class documentation
        """
        Plugin.setup(self, kwargs)
        self.c216 = self.input.get("c216", "id02/c216/0")
        self.cycle = self.input.get("cycle", 1)
        if "hdf5_filename" not in self.input:
            self.log_error("hdf5_filename not in input")
        self.hdf5_filename = self.input.get("hdf5_filename")
        self.entry = self.input.get("entry", "entry")
        self.instrument = self.input.get("instrument", "id02")

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

        # MultiCounterScaler
        self.mcs_grp = self.hdf5.require_group(posixpath.join(self.instrument, "MCS"))
        self.mcs_grp.attrs["NX_class"] = "NXcollection"
        self.mcs_grp["device"] = numpy.string_(self.c216)

        # Static metadata
        self.info_grp = self.hdf5.require_group(posixpath.join(self.instrument, "Information"))
        self.info_grp.attrs["NX_class"] = "NXcollection"
#        fields = ("MachineInfo", "OpticsInfo", "ProposalInfo", "StationInfo", "DetectorInfo", "ExperimentInfo") + \
#                 self.input.get("info", ())
        for field, value in self.input.get("info", {}).iteritems():
            self.info_grp[field] = numpy.string_(value)

        start_time = self.input.get("HMStartTime", get_isotime())

        # Factor
        HS32F = self.input.get("HS32F")
        if HS32F is not None:
            self.mcs_grp["HS32F"] = HS32F
        # Zero
        HS32Z = self.input.get("HS32Z")
        if HS32Z is not None:
            self.mcs_grp["HS32Z"] = HS32Z
        # Name
        HS32N = self.input.get("HS32N")
        if HS32N is not None:
            self.mcs_grp["HS32N"] = [str(i) for i in HS32N]
        # Mode
        HS32M = self.input.get("HS32M")
        if HS32M is not None:
            self.mcs_grp["HS32M"] = HS32M

        if HS32N and HS32Z and HS32F:
            self.mcs_grp.require_group("interpreted")
        self.group.parent["title"] = numpy.string_("id02.metadata")
        self.group.parent["program"] = numpy.string_("Dahu")
        self.group.parent["start_time"] = numpy.string_(start_time)

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
        # handle the case of multiple cycles: n*more frames, exposure time always the same
        tmp = numpy.outer(numpy.ones(cycles), tfg).ravel()
        live_time = tmp[1::2]
        self.tfg_grp["live_time"] = live_time
        self.tfg_grp["live_time"].attrs["interpretation"] = "scalar"
        self.tfg_grp["dead_time"] = tmp[0::2]
        self.tfg_grp["dead_time"].attrs["interpretation"] = "scalar"
        self.tfg_grp["delta_time"] = tmp.cumsum()[0::2]
        # raw scalers:
        raw_scalers = c216ds.ReadScalersForNLiveFrames([0, frames - 1])
        raw_scalers.shape = frames, -1
        counters = raw_scalers.shape[1]
        self.mcs_grp["HS32C"] = raw_scalers
        if "HSTime" in self.input:
            pin = int(self.input["HSTime"])
            if pin > counters:
                self.log_error("invalid pin number %s" % pin)
            self.mcs_grp["HSTime"] = pin
            self.mcs_grp["HSTime"].attrs["interpretation"] = "scalar"
            self.mcs_grp["HSTime"].attrs["counter"] = "1-based pin number"
            pin -= 1  # 1 based pin number
            time_counter = raw_scalers[:, pin]
            if "HS32F" in self.mcs_grp:
                factor = self.mcs_grp["HS32F"][pin]
            else:
                self.log_error("No factors provided for time measurement: defaulting to 1e-6", False)
                factor = 1e-6
            measured_time = time_counter * factor
            self.mcs_grp["meas_time"] = measured_time
            self.mcs_grp["meas_time"].attrs["interpretation"] = "scalar"
        else:
            self.log_error("No HSTime pin number, using TFG time")
            measured_time = tfg[1::2]

        if ("HS32F" in self.mcs_grp) and ("HS32Z" in self.mcs_grp):
            for I in ("HSI0", "HSI1"):
                if I in self.input:
                    pin = int(self.input[I])
                    if pin > counters:
                        self.log_error("invalid pin number %s" % pin)
                    self.mcs_grp[I] = pin
                    self.mcs_grp[I].attrs["interpretation"] = "scalar"
                    self.mcs_grp[I].attrs["counter"] = "1-based pin number"
                    pin -= 1  # 1 based pin number
                    counter = raw_scalers[:, pin]
                    factor = self.mcs_grp["HS32F"][pin]
                    zero = self.mcs_grp["HS32Z"][pin]
                    measured = (counter - measured_time * zero) * factor
                    self.mcs_grp[I[2:]] = measured
                    self.mcs_grp[I[2:]].attrs["interpretation"] = "scalar"
        else:
            self.log_error("Not factor/zero to calculate I0/I1", True)

        if "interpreted" in self.mcs_grp:
            modes = numpy.zeros(counters, dtype=numpy.int32)
            if "HS32M" in self.mcs_grp:
                raw_modes = numpy.array(self.mcs_grp["HS32M"])
                modes[:raw_modes.size] = raw_modes
            # else: mode=0
            values = numpy.zeros((frames, counters), dtype=numpy.float32)
            exptime = numpy.outer(tfg[1::2], numpy.ones(counters))
            zero = numpy.outer(numpy.ones(frames), numpy.array(self.mcs_grp["HS32Z"]))
            factor = numpy.outer(numpy.ones(frames), numpy.array(self.mcs_grp["HS32F"]))
            values_int = (raw_scalers - zero * exptime) * factor
            values_avg = (raw_scalers / exptime - zero) * factor
            mask = (numpy.outer(numpy.ones(frames), modes)).astype(bool)
            values[mask] = values_int[mask]
            nmask = numpy.logical_not(mask)
            values[nmask] = values_avg[nmask]
            self.mcs_grp["HS32V"] = values.astype(numpy.float32)
            self.mcs_grp["HS32V"].attrs["interpretation"] = "scalar"
            for i, name in enumerate(self.mcs_grp["HS32N"]):
                fullname = "interpreted/%s" % name
                self.mcs_grp[fullname] = values[:, i]
                self.mcs_grp[fullname].attrs["interpretation"] = "scalar"

    def teardown(self):
        if self.group:
            self.group.parent["end_time"] = numpy.string_(get_isotime())
        if self.hdf5:
            self.hdf5.close()
        Plugin.teardown(self)

################################################################################
# Image filtering plugin
################################################################################
@register
class Filter(Plugin):
    """
    This plugin does filtering of a set of images taken in HDF5 and outputs a single image file (in TIF or EDF)

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
        pyFAI.utils.averageImages(self.read_stack(), filter_=self.filter, cutoff=self.cutoff,
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
                        return numpy.array(ds)
                    
        
################################################################################
# Single Detector plugin
################################################################################
@register
class SingleDetector(Plugin):
    """
    This plugin does all processing needed for a single camera
    input = { "DetectorName": "rayonix",
              "image_file": "/nobackup/lid02gpu11/FRELON/test_laurent_saxs_0000.h5",
              #"entry": "entry_0000"
              #"hdf5": "/entry_0000/id02/data,
              "output_dir": "/nobackup/lid02gpu12",
              "PSize_1": 2.4e-05,
              "PSize_2" 2.4e-05,
              "BSize_1":1,
              "BSize_2":1,
              "Center_1":512,
              "Center_2":512,
              "DDummy":1,
              "SampleDistance":14.9522,
              "c216_filename": "/nobackup/lid02gpu11/metadata/test.h5"
              "WaveLength": 9.95058e-11,
              "Dummy":-10,
              "output_dir: "/nobackup/lid02gpu12/output",
#              "do_dark":false,
#              "do_flat":false,
              "npt2_azim": 360,
              "npt2_rad" : 500,
              "npt1_rad" : 1000,
              "to_save": ["raw", "dark", "flat", "dist", "norm", "azim", "ave"]
              }  
    """
    KEY_CONV = {"BSize": int,
                "Offset": int,
                "RasterOrientation": int,
                "Center": float,
                "DDummy": float,
                "Dummy": float,
                "PSize": float,
                "SampleDistance": float,
                "Wavelength": float,
                "Rot": float
                }

    KEYS = ("BSize_1", "BSize_2", "Center_1", "Center_2", "DDummy", "DetectorName",
            "Dummy", "Offset_1", "Offset_2", "PSize_1", "PSize_2",
            "Rot_1", "Rot_2", "Rot_3",
            "RasterOrientation", "SampleDistance", "SaxsDataVersion", "Title", "WaveLength")

    def __init__(self):
        Plugin.__init__(self)
        self.ai = None
        self.distotion_cor = None
        self.distotion_norm = None
        self.workers = {}
        self.output_ds = {} #output datasets
        self.dest = None    # output directory
        self.I1 = None      # beam stop diode values
        self.nframes = None
        self.to_save = ["raw", "ave"]  # by default only raw image and averaged one is saved
        self.input_nxs = None
        self.images_ds = None
        self.metadata = {}
        self.npt1_rad = None
        self.npt2_rad = None
        self.npt2_azim = None
        self.dark = None
        self.flat = None
        self.output_hdf5 = {}

    def setup(self, kwargs=None):
        """
        see class documentation
        """
        Plugin.setup(self, kwargs)
        if "output_dir" not in self.input:
            self.log_error("output_dir not in input")
        self.dest = os.path.abspath(self.input["output_dir"])
        if not os.path.isdir(self.dest):
            os.makedirs(self.dest)
        c216_filename = os.path.abspath(self.input.get("c216_filename", ""))
        if os.path.dirname(c216_filename) != self.dest:
            self.output_hdf5["metadata"] = os.path.join(self.dest, os.path.basename(c216_filename))
            t = threading.Thread(target=shutil.copy, name="copy metadata", args=(c216_filename, self.dest))
            t.start()
#            shutil.copy(c216_filename, self.dest)

        if "to_save" in self.input:
            self.to_save = self.input["to_save"][:]

        if "image_file" not in self.input:
            self.log_error("image_file not in input")
        self.image_file = self.input["image_file"]
        if not os.path.exists(self.image_file):
            self.log_error("image_file %s does not exist" % self.image_file)
        if "raw" in self.to_save:
            t = threading.Thread(target=shutil.copy, name="copy raw", args=(self.image_file, self.dest))
            t.start()
            self.to_save.remove("raw")
            self.output_hdf5["raw"] = os.path.join(self.dest, os.path.basename(self.image_file))
        self.hdf5_filename = self.input.get("hdf5_filename")
        self.entry = self.input.get("entry", "entry")
        self.instrument = self.input.get("instrument", "id02")
        self.I1 = self.load_I1(c216_filename)

    def process(self):
        self.metadata = self.parse_image_file()
        if self.I1 is None:
            self.I1 = numpy.ones(self.images_ds.shape[0], dtype=float)
        elif self.I1.size < self.images_ds.shape[0]:
            ones = numpy.ones(self.images_ds.shape[0], dtype=float)
            ones[:self.I1.size] = self.I1
            self.I1 = ones
        # update all metadata with the one provided by input
        for key, value in self.input.iteritems():
            if key in self.KEYS:
                self.metadata[key] = value
        forai = {}
        for key in ("BSize_1", "BSize_2", "Center_1", "Center_2",
                    "PSize_1", "PSize_2", "Rot_1", "Rot_2", "Rot_3",
                    "SampleDistance", "WaveLength"):
            forai[key] = self.metadata.get(key)
        self.ai = pyFAI.AzimuthalIntegrator()
        self.ai.setSPD(**forai)
        
        self.create_hdf5()
        self.process_images()

    def load_I1(self, mfile):
        """
        load the I1 data from a metadata HDF5 file

        /entry_0001/id02/MCS/I1

        @param mfile: metadata HDF5 file
        @return: array with I1
        """
        if "I1" in self.input:
            return numpy.array(self.input["I1"])
        nxs = pyFAI.io.Nexus(mfile, "r")
        for entry in nxs.get_entries():
            for instrument in nxs.get_class(entry, "NXinstrument"):
                if "MCS" in instrument:
                    mcs = instrument["MCS"]
                    if "I1" in mcs:
                        return numpy.array(mcs["I1"])
                
    def parse_image_file(self):
        """
        @return: dict with interpreted metadata
        """
        metadata = {}

        self.input_nxs = pyFAI.io.Nexus(self.image_file, "r")
        if "entry" in self.input:
            self.entry = self.input_nxs.get_entry(self.input["entry"])
        else:
            self.entry = self.input_nxs.get_entries()[0] #take the last entry
        instrument = self.input_nxs.get_class(self.entry, class_type="NXinstrument")
        if len(instrument) == 1:
            instrument = instrument[0]
        else:
            self.logg_error("Expected ONE instrument is expected in entry, got %s in %s %s" %
                            (len(instrument), self.image_file, self.entry))
        detector_grp = self.input_nxs.get_class(instrument, class_type="NXdetector")
        if len(detector_grp) == 1:
            detector_grp = detector_grp[0]
        elif len(detector_grp) == 0 and "detector" in instrument:
            detector_grp = instrument["detector"]
        else:
            self.logg_error("Expected ONE deteector is expected in experiment, got %s in %s %s %s" %
                            (len(detector_grp), self.input_nxs, self.image_file, instrument))
        self.images_ds = detector_grp.get("data")
        if "detector_information" in detector_grp:
            detector_information = detector_grp["detector_information"]
            if "name" in detector_information:
                metadata["DetectorName"] = str(detector_information["name"])
        #now read an interpret all static metadata.
        #metadata are on the collection side not instrument
        collections = self.input_nxs.get_class(self.entry, class_type="NXcollection")
        if len(collections) == 1:
            collection = collections[0]
        else:
            if len(collections) >= 1:
                collection = collections[0]
            else:
                self.logg_error("Expected ONE collections is expected in entry, got %s in %s %s" %
                            (len(collections), self.image_file, self.entry))

        detector_grps = self.input_nxs.get_class(collection, class_type="NXdetector")
        if len(detector_grps)==0 and "detector" in collection:
            detector_grp = collection["detector"]
        elif len(detector_grps) > 0:
            detector_grp = detector_grps[0]
        else:
            return {}
        if "parameters" in detector_grp:
            parameters_grp = detector_grp["parameters"]
        else:
            return {}

        for key, value in parameters_grp.iteritems():
            base = key.split("_")[0] 
            conv = self.KEY_CONV.get(base, str)
            metadata[key] = conv(value)
        return metadata
            
    def create_hdf5(self):
        """
        Create one HDF5 file per output
        Also initialize workers
        """
        in_shape = self.images_ds.shape
        for ext in self.to_save:
            outfile = os.path.join(self.dest, "%s_%s.h5" % (self.metadata["DetectorName"], ext))
            self.output_hdf5[ext] = outfile
            try:
                nxs = pyFAI.io.Nexus(outfile, "a")
            except IOError as error:
                self.log_error("invalid HDF5 file %s: remove and re-create!\n%s" % (outfile, error), False)
                os.unlink(outfile)
                nxs = pyFAI.io.Nexus(outfile)
            entry = nxs.new_entry("entry")
            subentry = nxs.new_class(entry, "pyFAI", class_type="NXsubentry")
            subentry["definition_local"] = numpy.string_("PyFAI")
            coll = nxs.new_class(subentry, "process_" + ext, class_type="NXcollection")
            metadata_grp = coll.require_group("metadata")
            for key, val in self.metadata.iteritems():
                if type(val) in [str, unicode]:
                    metadata_grp[key] = numpy.string_(val)
                else:
                    metadata_grp[key] = val
            if ext in ("raw", "dark", "flat", "cor", "norm"):
                shape = in_shape
            elif ext == "azim":
                if "npt2_rad" in self.input:
                    self.npt2_rad = int(self.input["npt2_rad"])
                else:
                    qmax = self.ai.qArray(in_shape[-2:]).max()
                    dqmin = self.ai.deltaQ(in_shape[-2:]).min() * 2.0
                    self.npt2_rad = int(qmax / dqmin)

                if "npt2_azim" in self.input:
                    self.npt2_azim = int(self.input["npt2_rad"])
                else:
                    chi = self.ai.chiArray(in_shape[-2:])
                    self.npt2_azim = int(numpy.degrees(chi.max() - chi.min()))
                shape = (in_shape[0], self.npt2_azim, self.npt2_rad)
                ai = pyFAI.AzimuthalIntegrator()
                ai.setPyFAI(**self.ai.getPyFAI())
                worker = pyFAI.worker.Worker(ai, in_shape[-2:], (self.npt2_azim, self.npt2_rad), "q_nm^-1")
                worker.output = "numpy"
                worker.method = "ocl_csr_gpu"
                self.workers[ext] = worker
            elif ext == "ave":
                if "npt1_rad" in self.input:
                    self.npt1_rad = int(self.input["npt1_rad"])
                else:
                    qmax = self.ai.qArray(in_shape[-2:]).max()
                    dqmin = self.ai.deltaQ(in_shape[-2:]).min() * 2.0
                    self.npt1_rad = int(qmax / dqmin)
                shape = (in_shape[0], self.npt1_rad)
                worker = pyFAI.worker.Worker(self.ai, in_shape[-2:], (1, self.npt1_rad), "q_nm^-1")
                worker.output = "numpy"
                worker.method = "ocl_csr_gpu"
                if self.flat:
                    worker.setFlatfieldFile(self.flat)
                if self.dark:
                    worker.setDarkcurrentFile(self.dark)
                self.workers[ext] = worker
            output_ds = coll.create_dataset("data", shape, "float32",
                                            chunks=(1,) + shape[1:],
                                            maxshape=(None,) + shape[1:])
            output_ds.attrs["NX_class"] = "NXdata"
            if ext == "ave":
                output_ds.attrs["interpretation"] = "spectrum"
            else:
                output_ds.attrs["interpretation"] = "image"
            self.output_ds[ext] = output_ds

    def process_images(self):
        """
        Here we process images....
        """
        for i in range(self.images_ds.shape[0]):
            data = self.images_ds[i]
            for meth in self.to_save:
                res = None
                ds = self.output_ds[meth]
                if meth == "dark":
                    pass  #  TODO
                elif meth == "flat":
                    pass  #  TODO
                elif meth == "cor":
                    res = self.distortion.correct(ds)
                elif meth == "norm":
                    res = self.distortion.correct(ds) / self.I1[i]
                elif meth == "azim":
                    res = self.workers[meth].process(data)
                    res /= self.I1[i]
                    if i == 0:
                        if "q" not in ds.parent:
                            ds.parent["q"] = self.workers[meth].radial
                            ds.parent["q"].attrs["unit"] = "q_nm^-1"
                        if "chi" not in ds.parent:
                            ds.parent["chi"] = self.workers[meth].azimuthal
                            ds.parent["chi"].attrs["unit"] = "deg"
                elif meth == "ave":
                    res = self.workers[meth].process(data)
                    if i == 0 and "q" not in ds.parent:
                        ds.parent["q"] = self.workers[meth].radial
                        ds.parent["q"].attrs["unit"] = "q_nm^-1"
                    res = res[:, 1]
                    res /= self.I1[i]
                ds[i] = res

    def teardown(self):
        if self.images_ds:
            self.images_ds.file.close()
        for ds in self.output_ds.values():
            ds.file.close()
        self.ai = None
        self.output["files"] = self.output_hdf5
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
