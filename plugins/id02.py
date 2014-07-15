#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement, print_function


import PyTango
import h5py
import logging
import numpy
import pyFAI.distortion
import sys, os, time, posixpath

from dahu.factory import register
from dahu.plugin import Plugin
from dahu.utils import get_isotime


__doc__ = """
Plugins for ID02:

* Distortion correction
* Metadata saving (C216)
* 
"""
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "10/072014"
__status__ = "development"
version = "0.2"

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

@register
class Metadata(Plugin):
    """
    Plugin in charge of retrieving all metadata for ID02 and storing them into a HDF5 file


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






    """       
ID02META_GENERAL["DetectorInfo"] = "VxH:detbox=14952.235x0.000x1.000,dettab=-62.000x-245.000"
ID02META_GENERAL["ExperimentInfo"] = 0
ID02META_GENERAL["HS32Depth"] = 32
ID02META_GENERAL["HS32F01"] = 1e-06
ID02META_GENERAL["HS32F02"] = 1
ID02META_GENERAL["HS32F03"] = 7763480
ID02META_GENERAL["HS32F04"] = 8176290
ID02META_GENERAL["HS32F05"] = 342239
ID02META_GENERAL["HS32F06"] = 967341
ID02META_GENERAL["HS32F07"] = 5541980
ID02META_GENERAL["HS32F08"] = 1739160
ID02META_GENERAL["HS32F09"] = 2753.61
ID02META_GENERAL["HS32F10"] = 1351920
ID02META_GENERAL["HS32F11"] = 140000000
ID02META_GENERAL["HS32F12"] = 16719500
ID02META_GENERAL["HS32F13"] = 1
ID02META_GENERAL["HS32F14"] = 0.000995868
ID02META_GENERAL["HS32F15"] = 0.000995868
ID02META_GENERAL["HS32F16"] = 1          
ID02META_GENERAL["HS32Len"] = 16
ID02META_GENERAL["HS32N01"] = "TIME"
ID02META_GENERAL["HS32N02"] = "AUX1"
ID02META_GENERAL["HS32N03"] = "PIN41"
ID02META_GENERAL["HS32N04"] = "PIN42"
ID02META_GENERAL["HS32N05"] = "PIN5"
ID02META_GENERAL["HS32N06"] = "PIN6"
ID02META_GENERAL["HS32N07"] = "PIN7"
ID02META_GENERAL["HS32N08"] = "PIN8"
ID02META_GENERAL["HS32N09"] = "PIN1"
ID02META_GENERAL["HS32N10"] = "PIN2"
ID02META_GENERAL["HS32N11"] = "PIN3"
ID02META_GENERAL["HS32N12"] = "PIN4"
ID02META_GENERAL["HS32N13"] = "AUX2"
ID02META_GENERAL["HS32N14"] = "THC1"
ID02META_GENERAL["HS32N15"] = "THC2"
ID02META_GENERAL["HS32N16"] = "PRESS"
ID02META_GENERAL["HS32Z01"] = 0
ID02META_GENERAL["HS32Z02"] = 0
ID02META_GENERAL["HS32Z03"] = 383.55
ID02META_GENERAL["HS32Z04"] = 126.4
ID02META_GENERAL["HS32Z05"] = 6582.1
ID02META_GENERAL["HS32Z06"] = 6973.6
ID02META_GENERAL["HS32Z07"] = 162.95
ID02META_GENERAL["HS32Z08"] = 0
ID02META_GENERAL["HS32Z09"] = 221.2
ID02META_GENERAL["HS32Z10"] = 207.8
ID02META_GENERAL["HS32Z11"] = 315.26
ID02META_GENERAL["HS32Z12"] = 145.11
ID02META_GENERAL["HS32Z13"] = 323.76
ID02META_GENERAL["HS32Z14"] = 170
ID02META_GENERAL["HS32Z15"] = 170
ID02META_GENERAL["HS32Z16"] = 228.4
ID02META_GENERAL["HSI0"] = 12
ID02META_GENERAL["HSI1"] = 7
ID02META_GENERAL["HSTime"] = 1
ID02META_GENERAL["MachineInfo"] = "Ie=183.27mA,u35u=100.000mm/0.001mm,u21m=100.000mm/0.000mm,u21d=100.000mm/-0.000mm"
ID02META_GENERAL["MirrorInfo"] = "rz=-3.600mrad,ty=0.455mm,ty1=2.075mm,ty2=-1.165mm,tz1=-0.030mm,tz2=-0.090mm,mir2rz=2.000mrad"
ID02META_GENERAL["OpticsInfo"] = "egy=12460.0eV,theta=9.132deg,hgt=11.7mm,tilt=4.440deg,tra=1.963mm"
ID02META_GENERAL["ProposalInfo"] = 0
ID02META_GENERAL["StationInfo"] = "ID02"
ID02META_GENERAL_NAME = "ID02META_GENERAL"



19.FRELON> syms -v ID02*STATIC**

ID02META_STATIC_NAME["frelon"] = "ID02META_STATIC_frelon"
ID02META_STATIC_NAME["saxs"] = "ID02META_STATIC_saxs"
ID02META_STATIC_frelon["BSize_1"] = 0
ID02META_STATIC_frelon["BSize_2"] = 0
ID02META_STATIC_frelon["Center_1"] = 0
ID02META_STATIC_frelon["Center_2"] = 0
ID02META_STATIC_frelon["DDummy"] = 0
ID02META_STATIC_frelon["DetectorName"] = 0
ID02META_STATIC_frelon["DetectorPosition"] = 14.9522
ID02META_STATIC_frelon["Dummy"] = 0
ID02META_STATIC_frelon["Offset_1"] = 0
ID02META_STATIC_frelon["Offset_2"] = 0
ID02META_STATIC_frelon["PSize_1"] = 0
ID02META_STATIC_frelon["PSize_2"] = 0
ID02META_STATIC_frelon["RasterOrientation"] = 0
ID02META_STATIC_frelon["SampleDistance"] = 14.9522
ID02META_STATIC_frelon["SaxsDataVersion"] = "saxs_data_version"
ID02META_STATIC_frelon["Title"] = 0
ID02META_STATIC_frelon["WaveLength"] = 9.95058e-11
 
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

@register
class BlaBla(Plugin):
    """
    This plugin does all processing needed
    input = { save = {"raw", "dark", "flat", "distortion", "normalization",  
    """
    def __init__(self):
        Plugin.__init__(self)

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
