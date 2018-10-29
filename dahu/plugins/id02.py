#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

"""Plugins for ID02:

* Distortion correction
* Metadata saving (C216)
* single detector processing
"""

from __future__ import with_statement, print_function, division

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "29/10/2018"
__status__ = "production"
version = "0.8"


import PyTango
import h5py
import logging
import numpy
import os
import posixpath
import pyFAI
import pyFAI.worker
import pyFAI.io
import pyFAI.utils
import fabio
import shutil
import sys
import time
import threading
import copy
import json
from collections import namedtuple
import dahu
from dahu.factory import register
from dahu.plugin import Plugin, plugin_from_function
from dahu.utils import get_isotime
from dahu.job import Job
from dahu.cache import DataCache
import numexpr

if sys.version_info < (3, 0):
    StringTypes = (str, unicode)
else:
    StringTypes = (str, bytes)

logger = logging.getLogger("dahu.id02")
# set loglevel at least at INFO
if logger.getEffectiveLevel() > logging.INFO:
    logger.setLevel(logging.INFO)
elif logger.getEffectiveLevel() < logging.root.level:
    logger.setLevel(logging.root.level)

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

CacheKey = namedtuple("CacheKey", ["ai", "mask", "shape"])


def preproc(**d):
    """Take a dict as input and forms a metadata structure as output
    @param: any dict
    """
    dd = d.copy()
    if "job_id" in dd:
        dd.pop("job_id")
    list_f = []
    list_n = []
    list_z = []
    HS32Len = dd.get("HS32Len", 16)
    HS32Depth = dd.get("HS32Depth", 32)
    HSI0Factor = dd.get("HSI0Factor", 1)
    HSI1Factor = dd.get("HSI1Factor", 1)
    # "0.005 s"
    if "ShutterOpeningTime" in dd:
        value = dd["ShutterOpeningTime"]
        if type(value) in StringTypes:
            ShutterOpeningTime = float(value.split()[0])
        else:
            ShutterOpeningTime = float(value)
    else:
        ShutterOpeningTime = 0
    if "ShutterClosingTime" in dd:
        value = dd["ShutterClosingTime"]
        if type(value) in StringTypes:
            ShutterClosingTime = float(value.split()[0])
        else:
            ShutterClosingTime = float(value)
    else:
        ShutterClosingTime = 0
    for ind in map(lambda x: 'HS32F' + '{0:02d}'.format(x), range(1, HS32Len + 1)):
            list_f.append(float(dd[ind]))
    for ind in map(lambda x: 'HS32N' + '{0:02d}'.format(x), range(1, HS32Len + 1)):
            list_n.append(dd[ind])
    for ind in map(lambda x: 'HS32Z' + '{0:02d}'.format(x), range(1, HS32Len + 1)):
            list_z.append(float(dd[ind]))

    info_dir = {}
    for info_ind in dd:
        if info_ind[0:2].find('HS') == 0:
            continue
        elif info_ind[0:2].find('HM') == 0:
            continue
        else:
            info_dir[info_ind] = dd[info_ind]

    final_dir = {"HS32Len": HS32Len,
                 "HS32Depth": HS32Depth,
                 "HSI0Factor": HSI0Factor,
                 "HSI1Factor": HSI1Factor,
                 "ShutterOpeningTime": ShutterOpeningTime,
                 "ShutterClosingTime": ShutterClosingTime,
                 'instrument': 'id02',
                 'c216': 'id02/c216/0',
                 'HS32F': list_f,
                 'HS32Z': list_z,
                 'HS32N': list_n,
                 'Info': info_dir}
    for key in ['HMStartEpoch', 'HMStartTime', "hdf5_filename", "entry", "HSTime", "HSI0", "HSI1"]:
        if key in dd:
            final_dir[key] = dd[key]
    return final_dir
plugin_from_function(preproc)


@register
class Metadata(Plugin):
    """Plugin in charge of retrieving all metadata for ID02 and storing them into a HDF5 file

    TODO: rewrite using Nexus class from pyFAI: shorter code

    NOTA: pin number are 1-based (I0, I1, time)

    Structure of the input data:
input = {
        "hdf5_filename":"/nobackup/lid02gpu11/metadata/test.h5",
        "entry": "entry",
        "instrument":"ESRF-ID02",
        "c216":"id02/c216/0",
        "HS32F": [1e-06, 1, 7763480, 8176290, 342239, 967341, 5541980, 1739160, 2753.61, 1351920, 140000000, 16719500, 1, 0.000995868, 0.000995868, 1],
        "HS32Z": [0, 0, 383.55, 126.4, 6582.1, 6973.6, 162.95, 0, 221.2, 207.8, 315.26, 145.11, 323.76, 170, 170, 228.4],
        "HS32N": ["TIME", "AUX1", "PIN41", "PIN42", "PIN5", "PIN6", "PIN7", "PIN8", "PIN1", "PIN2", "PIN3", "PIN4", "AUX2", "THC1", "THC2", "PRESS"],
        "HSI0": 12,
        "HSI1": 7,
        "HSTime": 1,
        "HMStartEpoch": 1405087717.12159,
        "HMStartTime": "2014-07-11T16:08:37.121591+0200",
        "Info": {"DetectorInfo":"VxH:detbox=14952.235x0.000x1.000,dettab=-62.000x-245.000",
                 "ExperimentInfo":"0",
                 "MachineInfo": "Ie=183.27mA,u35u=100.000mm/0.001mm,u21m=100.000mm/0.000mm,u21d=100.000mm/-0.000mm",
                 "MirrorInfo": "rz=-3.600mrad,ty=0.455mm,ty1=2.075mm,ty2=-1.165mm,tz1=-0.030mm,tz2=-0.090mm,mir2rz=2.000mrad",
                 "OpticsInfo": "egy=12460.0eV,theta=9.132deg,hgt=11.7mm,tilt=4.440deg,tra=1.963mm",
                 "ProposalInfo": 0,
                 "StationInfo": "ID02"
                 }
        }

"""
    TO_SKIP = ("entry", "hdf5_filename", "plugin_name")

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
        self.input2 = {}

    def setup(self, kwargs=None):
        """
        see class documentation
        """
        Plugin.setup(self, kwargs)
        if "HS32F10" in self.input:
            self.input2.update(preproc(**self.input))
        else:
            self.input2.update(self.input)
        # for debugging
        self.input["input2"] = self.input2

        self.c216 = self.input2.get("c216", "id02/c216/0")
        self.cycle = self.input2.get("cycle", 1)
        if "hdf5_filename" not in self.input2:
            self.log_error("hdf5_filename not in input")
        self.hdf5_filename = self.input2.get("hdf5_filename")
        self.entry = self.input2.get("entry", "entry")
        self.instrument = self.input2.get("instrument", "ESRF-ID02")

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
            self.log_warning("Unable to open %s: %s. Removing file and starting from scratch" % (self.hdf5_filename, error))
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
        self.info_grp = self.hdf5.require_group(posixpath.join(self.instrument, "parameters"))
        self.info_grp.attrs["NX_class"] = "NXcollection"

        for field, value in self.input2.get("Info", {}).items():
            if field not in self.TO_SKIP and not isinstance(value, dict):
                try:
                    value.encode("ascii")
                except UnicodeEncodeError:
                    self.log_warning("Unicode Error in field %s: %s, skipping" % (field, value))
                except AttributeError as err:
                    self.log_warning("Attribute Error %s \n in field %s: %s, forcing to string." % (err, field, value))
                    self.info_grp[field] = numpy.string_(value)
                else:
                    self.info_grp[field] = numpy.string_(value)

        start_time = self.input2.get("HMStartTime", get_isotime())

        # Factor
        HS32F = self.input2.get("HS32F")
        if HS32F is not None:
            self.mcs_grp["HS32F"] = HS32F
        # Zero
        HS32Z = self.input2.get("HS32Z")
        if HS32Z is not None:
            self.mcs_grp["HS32Z"] = HS32Z
        # Name
        HS32N = self.input2.get("HS32N")
        if HS32N is not None:
            self.mcs_grp["HS32N"] = numpy.array([str(i) for i in HS32N])
        # Mode
        HS32M = self.input2.get("HS32M")
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
        status = {"TFU_MODE": raw_status[0],
                  "TFU_FRAME_MODE": raw_status[1],
                  "START_MODE": raw_status[2],
                  "FRAMES": raw_status[14] // 2,
                  "CYCLES": raw_status[17],
                  "CHANNELS": raw_status[19],
                  "ACQ_BANK": raw_status[65],
                  "TFU_CONFIG": raw_status[:28],
                  "DIGI_IN": raw_status[28:44],
                  "ANA_IN": raw_status[44:60]}

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
        self.tfg_grp["delta_time"].attrs["interpretation"] = "scalar"
        for key in ["HMStartTime", "HMStartEpoch"]:
            if key in self.input2:
                if type(self.input2[key]) in StringTypes:
                    self.tfg_grp[key] = numpy.string_(self.input2[key])
                else:
                    self.tfg_grp[key] = self.input2[key]
                self.tfg_grp[key].attrs["interpretation"] = "scalar"

        # raw scalers:
        raw_scalers = c216ds.ReadScalersForNLiveFrames([0, frames - 1])
        raw_scalers.shape = frames, -1
        counters = raw_scalers.shape[1]
        self.mcs_grp["HS32C"] = raw_scalers
        if "HSTime" in self.input2:
            pin = int(self.input2["HSTime"])
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
                self.log_warning("No factors provided for time measurement: defaulting to 1e-6")
                factor = 1e-6
            measured_time = time_counter * factor
            self.mcs_grp["ExposureTime"] = measured_time
            self.mcs_grp["ExposureTime"].attrs["interpretation"] = "scalar"
        else:
            self.log_error("No HSTime pin number, using TFG time")
            measured_time = tfg[1::2]

        if ("HS32F" in self.mcs_grp) and ("HS32Z" in self.mcs_grp):
            #             if "interpreted" in self.mcs_grp:
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
            nmask = numpy.logical_not(mask)
            values[mask] = values_avg[mask]
            values[nmask] = values_int[nmask]
            self.mcs_grp["HS32V"] = values.astype(numpy.float32)
            self.mcs_grp["HS32V"].attrs["interpretation"] = "scalar"
            for i, name in enumerate(self.mcs_grp["HS32N"]):
                fullname = "interpreted/%s" % name
                self.mcs_grp[fullname] = values[:, i]
                self.mcs_grp[fullname].attrs["interpretation"] = "scalar"

            sot = self.input2.get("ShutterOpeningTime", 0.0)
            sct = self.input2.get("ShutterClosingTime", 0.0)
            for name, value in (("ShutterOpeningTime", sot),
                                ("ShutterClosingTime", sct)):
                        self.mcs_grp[name] = value
                        self.mcs_grp[name].attrs["interpretation"] = "scalar"
            correction_time = (measured_time - sot + sct) / (measured_time - sot)

            for I in ("HSI0", "HSI1"):
                if I in self.input2:
                    dest = "Intensity" + I[-1]
                    pin = int(self.input2[I])
                    if pin > counters:
                        self.log_error("invalid pin number %s" % pin)
                    self.mcs_grp[I] = pin
                    self.mcs_grp[I].attrs["interpretation"] = "scalar"
                    self.mcs_grp[I].attrs["counter"] = "1-based pin number"
                    pin -= 1  # 1 based pin number got 0 based.
                    counter = values[:, pin]
#                     factor = self.mcs_grp["HS32F"][pin]
#                     zero = self.mcs_grp["HS32Z"][pin]
#                     measured = (counter - measured_time * zero) * factor
                    I_factor = float(self.input2.get(I + "Factor", 1.0))
                    self.mcs_grp[I + "Factor"] = I_factor
                    self.mcs_grp[I + "Factor"].attrs["interpretation"] = "scalar"
                    measured = counter * I_factor
                    self.mcs_grp[dest] = measured
                    self.mcs_grp[dest].attrs["interpretation"] = "scalar"
                    self.mcs_grp[dest + "ShutCor"] = measured * correction_time
                    self.mcs_grp[dest + "ShutCor"].attrs["interpretation"] = "scalar"
        else:
            self.log_error("Not factor/zero to calculate I0/I1")

    def teardown(self):
        self.output["c216_filename"] = self.hdf5_filename
        if self.group:
            self.output["c216_path"] = self.group.name
            self.group.parent["end_time"] = numpy.string_(get_isotime())
        if self.hdf5:
            self.hdf5.close()
        Plugin.teardown(self)


################################################################################
# Single Detector plugin
################################################################################
@register
class SingleDetector(Plugin):
    """This plugin does all processing needed for a single camera

Minimalistic example:
              {
 "npt1_rad": 1000,
 "c216_filename": "/nobackup/lid02gpu11/metadata/test.h5",
 "npt2_rad": 500,
 "DetectorName": "rayonix",
 "npt2_azim": 360,
 "to_save": "raw sub ave azim dist flat ave_log",
 "output_dir": "/nobackup/lid02gpu12/output",
 "WaveLength": 9.95058e-11,
 "image_file": "/nobackup/lid02gpu11/FRELON/test_laurent_saxs_0000.h5",
 "dark_filename": "/nobackup/lid02gpu11/FRELON/test_laurent_dark_0000.h5",
}

All possible options are :
Mandatory options:
"image_file": HDF5 file containing the raw data
"output_dir": directory where final data are saved
"to_save": list of processing to be performed like: "raw sub flat dist norm azim ave"
"c216_filename": File containing the metadata (I1, timimgs ...)
"npt2_azim": number on bins in q direction for azim file
"npt2_rad": number on bins in azimuthal direction for azim file
"npt1_rad": number of bins in q direction for ave file

Optional parameters:
"DetectorName": "rayonix",
"PSize_1": 2.4e-05,
"PSize_2": 2.4e-05,
"BSize_1":1,
"BSize_2":1,
"Center_1":512,
"Center_2":512,
"DDummy":1,
"SampleDistance":14.9522,
"WaveLength": 9.95058e-11,
"Dummy":-10,
"metadata_job": int: wait for this job to finish and retrieve metadata from it.
"regrouping_mask_filename": file containing the mask to be used to integration
"dark_filename" HDF5 file with dark data in it
'dark_filter"  can be 'quantil' to average between given lower and upper quantils
"dark_filter_quantil" ask for the median id 0.5, min if 0 or max if 1
"dark_filter_quantil_lower" lower quantil for averaging
"dark_filter_quantil_upper" upper quantil for averaging
"distortion_filename" Spline file with the distortion correction
"flat_filename" flat field for intensity response correction
"metadata_job": number of the metadata job to wait for (unused)
"scaling_factor": float (default=1) by which all intensity will be multiplied 
"correct_solid_angle": True by default, set to 0/False to disable such correction 
"correct_I1": True by default, set to false to deactivate scaling by Exposure time / transmitted intensity
"unit": "q_nm^-1" can be changed to "log(q)_m" for log(q) units
"variance_formula": calculate the variance from a formula like '0.1*(data-dark)+1.0' "

Unused and automatically generated:
"plugin_name':"id02.singledetector',
"job_id': 29,

Defined but yet unused keywords:
"RasterOrientation": int,
"Rot": float

Possible values for to_save:
----------------------------

* sub: subtract dark signal
* flat: subtract dark then divide by flat
* solid: subtract dark then divide by flat and absolute solid angle
* dist: subtract dark then divide by flat and absolute solid angle, finally rebin on a regular grid
* norm: subtract dark then divide by flat and absolute solid angle, finally rebin on a regular grid and divide intensity by "I1"
* azim: subtract dark then divide by flat and absolute solid angle, finally rebin on a regular q/Chi grid and divide intensity by "I1"
* ave: subtract dark then divide by flat and absolute solid angle, finally rebin on a regular q grid and divide intensity by "I1"

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
    TIMEOUT = 10
    cache = DataCache(5)

    def __init__(self):
        Plugin.__init__(self)
        self.ai = None
        self.distotion_cor = None
        self.distotion_norm = None
        self.workers = {}
        self.output_ds = {}  # output datasets
        self.dest = None  # output directory
        self.I1 = None  # beam stop diode values
        self.t = None  # time of the start of the frame. Same shape as self.I1
        self.nframes = None
        self.to_save = ["raw", "ave"]  # by default only raw image and averaged one is saved
        self.input_nxs = None
        self.metadata_nxs = None
        self.images_ds = None
        self.metadata_plugin = None
        self.metadata = {}
        self.npt1_rad = None
        self.npt2_rad = None
        self.npt2_azim = None
        self.dark = None
        self.dark_filename = None
        self.flat_filename = None
        self.flat = None
        self.mask_filename = None
        self.distortion_filename = None
        self.output_hdf5 = {}
        self.dist = 1.0
        self.absolute_solid_angle = None
        self.in_shape = None
        self.scaling_factor = 1.0
        self.correct_solid_angle = True
        self.correct_I1 = True
        self.dummy = None
        self.delta_dummy = None
        self.unit = "q_nm^-1"
        self.polarization = None
        self.cachekey = None
        self.variance_formula = None
        self.variance_function = lambda data, dark: None

    def setup(self, kwargs=None):
        """
        see class documentation
        """
        Plugin.setup(self, kwargs)
        if "output_dir" not in self.input:
            self.log_error("output_dir not in input")
        self.dest = os.path.abspath(self.input["output_dir"])

        if "unit" in self.input:
            self.unit = self.input.get("unit")

        if "metadata_job" in self.input:
            job_id = int(self.input.get("metadata_job"))
            status = Job.synchronize_job(job_id, self.TIMEOUT)
            abort_time = time.time() + self.TIMEOUT
            while status == Job.STATE_UNINITIALIZED:
                # Wait for job to start
                time.sleep(1)
                status = Job.synchronize_job(job_id, self.TIMEOUT)
                if time.time() > abort_time:
                    self.log_error("Timeout while waiting metadata plugin to finish")
                    break
            if status == Job.STATE_SUCCESS:
                self.metadata_plugin = Job.getJobFromId(job_id)
            else:
                self.log_error("Metadata plugin ended in %s: aborting myself" % status)
        if not os.path.isdir(self.dest):
            os.makedirs(self.dest)
        c216_filename = os.path.abspath(self.input.get("c216_filename", ""))
        if (os.path.dirname(c216_filename) != self.dest) and (os.path.basename(c216_filename) not in os.listdir(self.dest)):
            self.output_hdf5["metadata"] = os.path.join(self.dest, os.path.basename(c216_filename))
            m = threading.Thread(target=shutil.copy, name="copy metadata", args=(c216_filename, self.dest))
            m.start()

        if "to_save" in self.input:
            to_save = self.input["to_save"][:]
            if type(to_save) in StringTypes:
                # fix a bug from spec ...
                self.to_save = [i.strip('[\\] ",') for i in to_save.split()]
                self.log_warning("processing planned: " + " ".join(self.to_save))
            else:
                self.to_save = to_save
        if "image_file" not in self.input:
            self.log_error("image_file not in input")
        self.image_file = self.input["image_file"]
        if not os.path.exists(self.image_file):
            if not self.image_file.startswith("/"):
                # prepend the dirname of the c216
                image_file = os.path.join(os.path.dirname(c216_filename), self.image_file)
                if os.path.exists(image_file):
                    self.image_file = image_file
                else:
                    self.log_error("image_file %s does not exist" % self.image_file)

        self.dark_filename = self.input.get("dark_filename")
        if "raw" in self.to_save:
            if os.path.dirname(self.image_file) != self.dest:
                t = threading.Thread(target=shutil.copy, name="copy raw", args=(self.image_file, self.dest))
                t.start()
            self.output_hdf5["raw"] = os.path.join(self.dest, os.path.basename(self.image_file))
            if type(self.dark_filename) in StringTypes and os.path.exists(self.dark_filename):
                if os.path.dirname(self.dark_filename) != self.dest:
                    d = threading.Thread(target=shutil.copy, name="copy dark", args=(self.dark_filename, self.dest))
                    d.start()
                self.output_hdf5["dark"] = os.path.join(self.dest, os.path.basename(self.dark_filename))
        self.scaling_factor = float(self.input.get("scaling_factor", 1.0))
        self.correct_solid_angle = bool(self.input.get("correct_solid_angle", True))
        self.correct_I1 = bool(self.input.get("correct_I1", True))
        self.I1, self.t = self.load_I1_t(c216_filename)

        # Variance formula: calculation of the function
        if "variance_formula" in self.input:
            self.variance_formula = self.input.get("variance_formula")
            if self.variance_formula:
                self.variance_function = numexpr.NumExpr(self.variance_formula,
                                                         [("data", float),
                                                          ("dark", float)])

    def process(self):
        self.metadata = self.parse_image_file()
        self.mask_filename = self.input.get("regrouping_mask_filename")
        shape = self.in_shape[-2:]
        if self.I1 is None:
            self.I1 = numpy.ones(shape, dtype=float)
        elif self.I1.size < self.in_shape[0]:
            ones = numpy.ones(self.in_shape[0], dtype=float)
            ones[:self.I1.size] = self.I1
            self.I1 = ones
        # update all metadata with the one provided by input
        for key, value in self.input.items():
            if key in self.KEYS:
                self.metadata[key] = value
        forai = {}
        for key in ("BSize_1", "BSize_2", "Center_1", "Center_2",
                    "PSize_1", "PSize_2", "Rot_1", "Rot_2", "Rot_3",
                    "SampleDistance", "WaveLength"):
            forai[key] = self.metadata.get(key)
        self.dist = self.metadata.get("SampleDistance")
        self.dummy = self.metadata.get("Dummy", self.dummy)
        self.delta_dummy = self.metadata.get("DDummy", self.delta_dummy)
        # read detector distortion distortion_filename
        self.ai = pyFAI.AzimuthalIntegrator()
        self.ai.setSPD(**forai)

        self.distortion_filename = self.input.get("distortion_filename") or None
        if type(self.distortion_filename) in StringTypes:
            detector = pyFAI.detectors.Detector(splineFile=self.distortion_filename)
            if tuple(detector.shape) != shape:
                self.log_warning("Binning needed for spline ? detector claims %s but image size is %s" % (detector.shape, shape))
                detector.guess_binning(shape)
            self.ai.detector = detector

        self.log_warning("AI:%s" % self.ai)

        self.cachekey = CacheKey(str(self.ai), self.mask_filename, shape)

        if self.cachekey in self.cache:
            self.ai._cached_array.update(self.cache.get(self.cachekey))
        else:
            # Initialize geometry:
            self.ai.qArray(shape)
            self.ai.chiArray(shape)
            self.ai.deltaQ(shape)
            self.ai.deltaChi(shape)
            self.ai.solidAngleArray(shape)
            self.cache.get(self.cachekey, {}).update(self.ai._cached_array)

        if self.input.get("do_polarization"):
            self.polarization = self.ai.polarization(factor=self.input.get("polarization_factor"),
                                                     axis_offset=self.input.get("polarization_axis_offset", 0))

        # Read and Process dark
        if type(self.dark_filename) in StringTypes and os.path.exists(self.dark_filename):
            dark = self.read_data(self.dark_filename)
            if dark.ndim == 3:
                method = self.input.get("dark_filter")
                if method.startswith("quantil"):
                    lower = self.input.get("dark_filter_quantil_lower", 0)
                    upper = self.input.get("dark_filter_quantil_upper", 1)
                    self.dark = pyFAI.average.average_dark(dark, center_method=method, quantiles=(lower, upper))
                else:
                    if method is None:
                        method = "median"
                    self.dark = pyFAI.average.average_dark(dark, center_method=method)
            else:
                self.dark = dark
        elif type(self.dark_filename) in (int, float):
            self.dark = float(self.dark_filename)
        self.ai.set_darkcurrent(self.dark)

        # Read and Process Flat

        self.flat_filename = self.input.get("flat_filename")
        if type(self.flat_filename) in StringTypes and os.path.exists(self.flat_filename):
            if self.flat_filename.endswith(".h5") or self.flat_filename.endswith(".nxs") or self.flat_filename.endswith(".hdf5"):
                flat = self.read_data(self.flat_filename)
            else:
                flat_fabio = fabio.open(self.flat_filename)
                flat = flat_fabio.data
                dummy = flat_fabio.header.get("Dummy")
                try:
                    dummy = float(dummy)
                except:
                    self.log_error("Dummy value in mask is unconsitent %s" % dummy)
                    dummy = None

            if flat.ndim == 3:
                self.flat = pyFAI.average.average_dark(flat, center_method="median")
            else:
                self.flat = flat
            if (self.flat is not None) and (self.flat.shape != shape):
                binning = [j / i for i, j in zip(shape, self.flat.shape)]
                if tuple(binning) != (1, 1):
                    self.log_warning("Binning for flat is %s" % binning)
                    if max(binning) > 1:
                        binning = [int(i) for i in binning]
                        self.flat = pyFAI.utils.binning(self.flat, binsize=binning, norm=False)
                    else:
                        binning = [i // j for i, j in zip(shape, self.flat.shape)]
                        self.flat = pyFAI.utils.unBinning(self.flat, binsize=binning, norm=False)
            self.ai.set_flatfield(self.flat)

        # Read and Process mask
        if type(self.mask_filename) in StringTypes and os.path.exists(self.mask_filename):
            try:
                mask_fabio = fabio.open(self.mask_filename)
            except:
                mask = self.read_data(self.mask_filename) != 0
            else:  # this is very ID02 specific !!!!
                dummy = mask_fabio.header.get("Dummy")
                try:
                    dummy = float(dummy)
                except:
                    self.log_error("Dummy value in mask is unconsitent %s" % dummy)
                    dummy = None
                ddummy = mask_fabio.header.get("DDummy")
                try:
                    ddummy = float(ddummy)
                except:
                    self.log_error("DDummy value in mask is unconsitent %s" % ddummy)
                    ddummy = 0
                if ddummy:
                    mask = abs(mask_fabio.data - dummy) <= ddummy
                else:
                    mask = (mask_fabio.data == dummy)
                self.log_warning("found %s pixel masked out" % (mask.sum()))
                self.dummy = dummy
                self.delta_dummy = ddummy
            if mask.ndim == 3:
                mask = pyFAI.average.average_dark(mask, center_method="median")
            if (mask is not None) and (mask.shape != shape):
                binning = [j / i for i, j in zip(shape, mask.shape)]
                if tuple(binning) != (1, 1):
                    self.log_warning("Binning for mask is %s" % binning)
                    if max(binning) > 1:
                        binning = [int(i) for i in binning]
                        mask = pyFAI.utils.binning(mask, binsize=binning, norm=True) > 0
                    else:
                        binning = [i // j for i, j in zip(shape, mask.shape)]
                        mask = pyFAI.utils.unBinning(mask, binsize=binning, norm=False) > 0
            # nota: this is assigned to the detector !
            self.ai.mask = mask

        # bug specific to ID02, dummy=0 means no dummy !
        if self.dummy == 0:
            self.dummy = None

        self.create_hdf5()
        self.process_images()

    def load_I1_t(self, mfile, correct_shutter_closing_time=True):
        """
        load the I1 data and timstamp for frame start from a metadata HDF5 file


        /entry_0001/id02/MCS/I1

        TODO: handle correction or not for shutter opening/closing time

        @param mfile: metadata HDF5 file
        @param correct_shutter_closing_time: set to true for integrating detector (CCD) and false for counting detector (Pilatus)
        @return: 2-tuple of array with I1 and t
        """
        if ("I1" in self.input):
            return numpy.array(self.input["I1"]), None
        self.metadata_nxs = pyFAI.io.Nexus(mfile, "r")
        I1 = t = None
        if correct_shutter_closing_time:
            key = "Intensity1ShutCor"
        else:
            key = "Intensity1"
        for entry in self.metadata_nxs.get_entries():
            for instrument in self.metadata_nxs.get_class(entry, "NXinstrument"):
                if "MCS" in instrument:
                    mcs = instrument["MCS"]
                    if key in mcs:
                        I1 = numpy.array(mcs[key])
                if "TFG" in instrument:
                    tfg = instrument["TFG"]
                    if "delta_time" in tfg:
                        t = numpy.array(tfg["delta_time"])
                if (t is None) or (I1 is None):
                    I1 = t = None
                else:
                    break
            if (t is not None) and (I1 is not None):
                return I1, t
        return I1, t

    def parse_image_file(self):
        """
        @return: dict with interpreted metadata
        """
        metadata = {}

        self.input_nxs = pyFAI.io.Nexus(self.image_file, "r")
        if "entry" in self.input:
            self.entry = self.input_nxs.get_entry(self.input["entry"])
        else:
            self.entry = self.input_nxs.get_entries()[0]  # take the last entry
        instrument = self.input_nxs.get_class(self.entry, class_type="NXinstrument")
        if len(instrument) == 1:
            instrument = instrument[0]
        else:
            self.log_error("Expected ONE instrument is expected in entry, got %s in %s %s" %
                           (len(instrument), self.image_file, self.entry))
        detector_grp = self.input_nxs.get_class(instrument, class_type="NXdetector")
        if len(detector_grp) == 1:
            detector_grp = detector_grp[0]
        elif len(detector_grp) == 0 and "detector" in instrument:
            detector_grp = instrument["detector"]
        else:
            self.log_error("Expected ONE deteector is expected in experiment, got %s in %s %s %s" %
                           (len(detector_grp), self.input_nxs, self.image_file, instrument))
        self.images_ds = detector_grp.get("data")
        if isinstance(self.images_ds, h5py.Group):
            if self.images_ds.attrs.get("NX_class") == "NXdata":
                # this is an NXdata not a dataset: use the @signal
                self.images_ds = self.images_ds.get(self.images_ds.attrs.get("signal"))
        self.in_shape = self.images_ds.shape
        if "detector_information" in detector_grp:
            detector_information = detector_grp["detector_information"]
            if "name" in detector_information:
                metadata["DetectorName"] = str(detector_information["name"])
        # now read an interpret all static metadata.
        # metadata are on the collection side not instrument
        collections = self.input_nxs.get_class(self.entry, class_type="NXcollection")
        if len(collections) == 1:
            collection = collections[0]
        else:
            if len(collections) >= 1:
                collection = collections[0]
            else:
                self.log_error("Expected ONE collections is expected in entry, got %s in %s %s" %
                               (len(collections), self.image_file, self.entry))

        detector_grps = self.input_nxs.get_class(collection, class_type="NXdetector")
        if (len(detector_grps) == 0) and ("detector" in collection):
            detector_grp = collection["detector"]
        elif len(detector_grps) > 0:
            detector_grp = detector_grps[0]
        else:
            return {}
        if "parameters" in detector_grp:
            parameters_grp = detector_grp["parameters"]
        else:
            return {}

        for key, value in parameters_grp.items():
            base = key.split("_")[0]
            conv = self.KEY_CONV.get(base, str)
            metadata[key] = conv(value.value)
        return metadata

    def create_hdf5(self):
        """
        Create one HDF5 file per output
        Also initialize all workers
        """
        basename = os.path.splitext(os.path.basename(self.image_file))[0]
        if basename.endswith("_raw"):
            basename = basename[:-4]
        json_config = json.dumps(self.input)
        isotime = numpy.string_(get_isotime())
        detector_grp = self.input_nxs.find_detector(all=True)
        detector_name = "undefined"
        for grp in detector_grp:
            if "detector_information/name" in grp:
                detector_name = grp["detector_information/name"].value
        md_entry = self.metadata_nxs.get_entries()[0]
        instruments = self.metadata_nxs.get_class(md_entry, "NXinstrument")
        if instruments:
            collections = self.metadata_nxs.get_class(instruments[0], "NXcollection")
            to_copy = collections + detector_grp
        else:
            to_copy = detector_grp

        for ext in self.to_save:
            if ext == "raw":
                continue
            outfile = os.path.join(self.dest, "%s_%s.h5" % (basename, ext))
            self.output_hdf5[ext] = outfile
            try:
                nxs = pyFAI.io.Nexus(outfile, "a")
            except IOError as error:
                self.log_warning("invalid HDF5 file %s: remove and re-create!\n%s" % (outfile, error))
                os.unlink(outfile)
                nxs = pyFAI.io.Nexus(outfile)
            entry = nxs.new_entry("entry", program_name="dahu", title=self.image_file + ":" + self.images_ds.name)

            entry["program_name"].attrs["version"] = dahu.version
            entry["plugin_name"] = numpy.string_(".".join((os.path.splitext(os.path.basename(__file__))[0], self.__class__.__name__)))
            entry["plugin_name"].attrs["version"] = version
            entry["input"] = numpy.string_(json_config)
            entry["input"].attrs["format"] = 'json'

            entry["detector_name"] = numpy.string_(detector_name)

            subentry = nxs.new_class(entry, "PyFAI", class_type="NXprocess")
            subentry["program"] = numpy.string_("PyFAI")
            subentry["version"] = numpy.string_(pyFAI.version)
            subentry["date"] = isotime
            subentry["processing_type"] = numpy.string_(ext)
            coll = nxs.new_class(subentry, "process_" + ext, class_type="NXdata")
            metadata_grp = coll.require_group("parameters")

            for key, val in self.metadata.iteritems():
                if type(val) in [str, unicode]:
                    metadata_grp[key] = numpy.string_(val)
                else:
                    metadata_grp[key] = val

            # copy metadata from other files:
            for grp in to_copy:
                grp_name = posixpath.split(grp.name)[-1]
                if grp_name not in coll:
                    toplevel = coll.require_group(grp_name)
                    for k, v in grp.attrs.items():
                        toplevel.attrs[k] = v
                else:
                    toplevel = coll[grp_name]

                def grpdeepcopy(name, obj):
                    nxs.deep_copy(name, obj, toplevel=toplevel, excluded=["data"])

                grp.visititems(grpdeepcopy)

            shape = self.in_shape[:]
            if self.npt1_rad is None and "npt1_rad" in self.input:
                self.npt1_rad = int(self.input["npt1_rad"])
            else:
                qmax = self.ai.qArray(self.in_shape[-2:]).max()
                dqmin = self.ai.deltaQ(self.in_shape[-2:]).min() * 2.0
                self.npt1_rad = int(qmax / dqmin)

            if ext == "azim":
                if "npt2_rad" in self.input:
                    self.npt2_rad = int(self.input["npt2_rad"])
                else:
                    qmax = self.ai.qArray(self.in_shape[-2:]).max()
                    dqmin = self.ai.deltaQ(self.in_shape[-2:]).min() * 2.0
                    self.npt2_rad = int(qmax / dqmin)

                if "npt2_azim" in self.input:
                    self.npt2_azim = int(self.input["npt2_azim"])
                else:
                    chi = self.ai.chiArray(self.in_shape[-2:])
                    self.npt2_azim = int(numpy.degrees(chi.max() - chi.min()))
                shape = (self.in_shape[0], self.npt2_azim, self.npt2_rad)

                ai = self.ai.__copy__()
                worker = pyFAI.worker.Worker(ai, self.in_shape[-2:], (self.npt2_azim, self.npt2_rad), self.unit)
                if self.flat is not None:
                    worker.ai.set_flatfield(self.flat)
                if self.dark is not None:
                    worker.ai.set_darkcurrent(self.dark)
                worker.output = "numpy"
                if self.in_shape[0] < 5:
                    worker.method = "splitbbox"
                else:
                    worker.method = "ocl_csr_gpu"
                if self.correct_solid_angle:
                    worker.set_normalization_factor(self.ai.pixel1 * self.ai.pixel2 / self.ai.dist / self.ai.dist)
                else:
                    worker.set_normalization_factor(1.0)
                    worker.correct_solid_angle = self.correct_solid_angle
                self.log_warning("Normalization factor: %s" % worker.normalization_factor)

                worker.dummy = self.dummy
                worker.delta_dummy = self.delta_dummy
                if self.input.get("do_polarization"):
                    worker.polarization_factor = self.input.get("polarization_factor")

                self.workers[ext] = worker
            elif ext.startswith("ave"):
                if "_" in ext:
                    unit = ext.split("_", 1)[1]
                    npt1_rad = self.input.get("npt1_rad_" + unit, self.npt1_rad)
                    ai = self.ai.__copy__()
                else:
                    unit = self.unit
                    npt1_rad = self.npt1_rad
                    ai = self.ai
                shape = (self.in_shape[0], npt1_rad)
                worker = pyFAI.worker.Worker(ai, self.in_shape[-2:], (1, npt1_rad), unit=unit)
                worker.output = "numpy"
                if self.in_shape[0] < 5:
                    worker.method = "splitbbox"
                else:
                    worker.method = "ocl_csr_gpu"
                if self.correct_solid_angle:
                    worker.set_normalization_factor(self.ai.pixel1 * self.ai.pixel2 / self.ai.dist / self.ai.dist)
                else:
                    worker.set_normalization_factor(1.0)
                    worker.correct_solid_angle = self.correct_solid_angle
                worker.dummy = self.dummy
                worker.delta_dummy = self.delta_dummy
                if self.input.get("do_polarization"):
                    worker.polarization_factor = True
                self.workers[ext] = worker
            elif ext == "sub":
                worker = pyFAI.worker.PixelwiseWorker(dark=self.dark,
                                                      dummy=self.dummy, delta_dummy=self.delta_dummy,
                                                      )
                self.workers[ext] = worker
            elif ext == "flat":
                worker = pyFAI.worker.PixelwiseWorker(dark=self.dark, flat=self.flat,
                                                      dummy=self.dummy, delta_dummy=self.delta_dummy,
                                                      )
                self.workers[ext] = worker
            elif ext == "solid":
                worker = pyFAI.worker.PixelwiseWorker(dark=self.dark, flat=self.flat, solidangle=self.get_solid_angle(),
                                                      dummy=self.dummy, delta_dummy=self.delta_dummy, polarization=self.polarization
                                                      )
                self.workers[ext] = worker
            elif ext == "dist":
                worker = pyFAI.worker.DistortionWorker(dark=self.dark, flat=self.flat, solidangle=self.get_solid_angle(),
                                                       dummy=self.dummy, delta_dummy=self.delta_dummy, polarization=self.polarization,
                                                       detector=self.ai.detector)
                self.workers[ext] = worker
            elif ext == "norm":
                worker = pyFAI.worker.DistortionWorker(dark=self.dark, flat=self.flat, solidangle=self.get_solid_angle(),
                                                       dummy=self.dummy, delta_dummy=self.delta_dummy, polarization=self.polarization,
                                                       detector=self.ai.detector)
                self.workers[ext] = worker
            else:
                self.log_warning("unknown treatment %s" % ext)

            # TODO: manage compression here

            output_ds = coll.create_dataset("data", shape, "float32",
                                            chunks=(1,) + shape[1:],
                                            maxshape=(None,) + shape[1:])
            coll.attrs["signal"] = "data"
            output_ds.attrs["signal"] = "1"
            if self.variance_formula is not None:
                error_ds = coll.create_dataset("errors", shape, "float32",
                                               chunks=(1,) + shape[1:],
                                               maxshape=(None,) + shape[1:])
                coll.attrs["uncertainties"] = "errors"
                self.output_ds[ext + "_err"] = error_ds
            if self.t is not None:
                coll["t"] = self.t
                coll["t"].attrs["axis"] = "1"
                coll["t"].attrs["interpretation"] = "scalar"
                coll["t"].attrs["unit"] = "s"

            if ext == "azim":
                coll.attrs["axes"] = ["t", "chi", "q"]
                output_ds.attrs["interpretation"] = "image"
            elif ext == "ave":
                coll.attrs["axes"] = ["t", "q"]
                output_ds.attrs["interpretation"] = "spectrum"
            elif ext in ("sub", "flat", "solid", "dist"):
                coll.attrs["axes"] = ["t", ".", "."]
                output_ds.attrs["interpretation"] = "image"
            else:
                output_ds.attrs["interpretation"] = "image"
            self.output_ds[ext] = output_ds

    def process_images(self):
        """
        Here we process images....
        """
        if self.correct_I1:
            I1s = self.I1
        else:
            I1s = numpy.ones_like(self.I1)
        for i, I1 in enumerate(I1s):
            data = self.images_ds[i]
            if self.variance_formula:
                variance = self.variance_function(data, self.dark)
            else:
                variance = None
            I1_corrected = I1 / self.scaling_factor
            for meth in self.to_save:
                if meth in ["raw", "dark"]:
                    continue
                res = err = None

                ds = self.output_ds[meth]
                if variance is not None:
                    err_ds = self.output_ds[meth + "_err"]

                if meth in ("sub", "flat", "dist", "cor"):
                    res = self.workers[meth].process(data, variance=variance)
                    if variance is not None:
                        res, err = res
                elif meth == "norm":
                    res = self.workers[meth].process(data, variance=variance,
                                                     normalization_factor=I1_corrected)
                    if variance is not None:
                        res, err = res

                elif meth == "azim":
                    res = self.workers[meth].process(data, variance=variance,
                                                     normalization_factor=I1_corrected)
                    if variance is not None:
                        res, err = res

                    if i == 0:
                        if "q" not in ds.parent:
                            ds.parent["q"] = numpy.ascontiguousarray(self.workers[meth].radial, dtype=numpy.float32)
                            ds.parent["q"].attrs["unit"] = self.unit
                            ds.parent["q"].attrs["axis"] = "3"
                            ds.parent["q"].attrs["interpretation"] = "scalar"
                        if "chi" not in ds.parent:
                            ds.parent["chi"] = numpy.ascontiguousarray(self.workers[meth].azimuthal, dtype=numpy.float32)
                            ds.parent["chi"].attrs["unit"] = "deg"
                            ds.parent["chi"].attrs["axis"] = "2"
                            ds.parent["chi"].attrs["interpretation"] = "scalar"
                elif meth.startswith("ave"):
                    res = self.workers[meth].process(data, variance=variance,
                                                     normalization_factor=I1_corrected)
                    if variance is not None:
                        res, err = res

                    # TODO: add other units
                    if i == 0 and "q" not in ds.parent:
                        if "log(1+q.nm)" in meth:
                            q = numpy.exp(self.workers[meth].radial) - 1.0
                        elif "log(1+q.A)" in meth:
                            q = (numpy.exp(self.workers[meth].radial) - 1.0) * 0.1
                        elif "log(q.nm)" in meth:
                            q = numpy.exp(self.workers[meth].radial)
                        elif "log10(q.m)" in meth:
                            q = 10 ** (self.workers[meth].radial) * 1e-9
                        else:
                            q = self.workers[meth].radial
                        ds.parent["q"] = numpy.ascontiguousarray(q, dtype=numpy.float32)
                        ds.parent["q"].attrs["unit"] = self.unit
                        ds.parent["q"].attrs["axis"] = "2"
                        ds.parent["q"].attrs["interpretation"] = "scalar"
                    res = res[:, 1]
                else:
                    self.log_warning("Unknown/supported method ... %s, skipping" % (meth))
                    continue
                ds[i] = res
                if err is not None:
                    err_ds[i] = err

    def read_data(self, filename):
        """read dark data from a file

        @param filename: HDF5 file containing dark frames
        @return: numpy array with dark
        """
        with pyFAI.io.Nexus(filename, "r") as nxs:
            for entry in nxs.get_entries():
                for instrument in nxs.get_class(entry, "NXinstrument"):
                    for detector in nxs.get_class(instrument, "NXdetector"):
                        if "data" in detector:
                            return numpy.array(detector["data"])

    def get_solid_angle(self):
        """ calculate the solid angle if needed and return it
        """
        if (self.absolute_solid_angle is None) and self.correct_solid_angle:
            self.absolute_solid_angle = self.ai.solidAngleArray(self.in_shape[-2:], absolute=True)
        return self.absolute_solid_angle

    def teardown(self):
        """Method for cleaning up stuff
        """
        # Finally update the cache:
        to_cache = {}
        for meth in self.to_save:
            worker = self.workers[meth]
            if "ai" in dir(worker):
                to_cache.update(worker.ai._cached_array)

        self.cache.get(self.cachekey, {}).update(to_cache)

        if self.images_ds:
            self.images_ds.file.close()
        for ds in self.output_ds.values():
            ds.file.close()
        self.ai = None
        self.polarization = None
        self.output["files"] = self.output_hdf5
        Plugin.teardown(self)
