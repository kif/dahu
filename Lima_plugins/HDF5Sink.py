#!/usr/bin/env python
# coding: utf8
from __future__ import with_statement, print_function
"""
LImA ProcessLib example of HDF5 writer 

This depends on PyFAI and on h5py

"""
__author__ = "Jérôme Kieffer"
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "GPLv3+"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "12/12/2013"
__status__ = "beta"
__docformat__ = 'restructuredtext'

import os, threading, logging, posixpath, time, types

import numpy
from Lima import Core
import h5py


import PyTango
import numpy
from Utils import BasePostProcess


class PrepareAcqCallback(Core.SoftCallback):
    """
    Class managing the connection from a 
    Lima.Core.CtControl.prepareAcq() to the configuration of the various tasks
    
    """

    CONFIG_ITEMS = {"dimX": None,
                    "dimY":None,
                    "binX": None,
                    "binY":None,
                    "directory": None,
                    "prefix": None,
                    "start_index": None,
                    "number_of_frames":None,
                    "exposure_time": None,
                     }
    LIMA_DTYPE = {Core.Bpp10:   "uint16",
                  Core.Bpp10S:   "int16",
                  Core.Bpp12:   "uint16",
                  Core.Bpp12S:   "int16",
                  Core.Bpp14:   "uint16",
                  Core.Bpp14S:   "int16",
                  Core.Bpp16:   "uint16",
                  Core.Bpp16S:   "int16",
                  Core.Bpp32:   "uint32",
                  Core.Bpp32S:   "int32",
                  Core.Bpp8:    "uint8",
                  Core.Bpp8S:   "int8"
                  }
    LIMA_ROTATION = {Core.Rotation_0: "no_rot",
                     Core.Rotation_90: "rot_90_cw",
                     Core.Rotation_270: "rot_90_ccw",
                     Core.Rotation_180: "rot_180",
                     }
    def __init__(self, control, task=None):
        """       
        @param control: Lima.Core.CtControl instance
        @param task: The task one wants to parametrize at startup. Can be a  Core.Processlib.LinkTask or a Core.Processlib.SinkTask
        """
        Core.SoftCallback.__init__(self)
        self._control = control
        self._task = task

    def prepare(self):
        """
        Called with prepareAcq()
        """

        im = self._control.image()
        imdim = im.getImageDim().getSize()
        x = imdim.getWidth()
        y = imdim.getHeight()
        bin = im.getBin()
        binX = bin.getX()
        binY = bin.getY()
        flip = im.getFlip()
        roi = im.getRoi()
        lima_cfg = {"dimX":x,
                    "dimY":y,
                    "binX":binX,
                    "binY":binY,
                    "flipX":flip.x,
                    "flipY":flip.y,
                    "rotation":self.LIMA_ROTATION[im.getRotation()],
                    "mode": im.getMode(),
                    "dtype": self.LIMA_DTYPE[im.getImageType()]}
        if roi.isActive():
            lima_cfg["OffsetX"] = roi.getTopLeft().x
            lima_cfg["OffsetY"] = roi.getTopLeft().y
        saving = self._control.saving()
        sav_parms = saving.getParameters()
        lima_cfg["directory"] = sav_parms.directory
        lima_cfg["prefix"] = sav_parms.prefix
        lima_cfg["start_index"] = sav_parms.nextNumber
        lima_cfg["indexFormat"] = sav_parms.indexFormat
        # number of images ...
        acq = self._control.acquisition()
        lima_cfg["number_of_frames"] = acq.getAcqNbFrames()
        lima_cfg["exposure_time"] = acq.getAcqExpoTime()

        if self._task._writer:
            self._task._writer.init(lima_cfg=lima_cfg)
            self._task._writer.flush()


class HDF5Sink(Core.Processlib.SinkTaskBase):
    """
    This is a ProcessLib task which is a sink: 
    it saves the image into a HDF5 stack.
    
    """

    def __init__(self, writer):
        Core.Processlib.SinkTaskBase.__init__(self)
        self._writer = writer

    def __repr__(self):
        """
        pretty print of myself
        """
        lstout = [ "HDF5Sink instance", "Writer:", self._writer.__repr__()]
        return os.linesep.join(lstout)


    def process(self, data) :
        """
        Callback function       
        Called for every frame in a different C++ thread.
        """
        self._writer.write(data.buffer, data.frameNumber)
        # horrible, only way to close the h5 file, but no garanty that all the thread
        # did finish to save !!
        if self._writer.isLastFrame(data.frameNumber):
            self._writer.close()


class HDF5Writer(object):
    """
    Class allowing to write HDF5 Files.    
    """
    CONFIG_ITEMS = {"filename": None,
                    "dirname":None,
                    "extension": ".h5",
                    "subdir":None,
                    "hpath": "entry_",
                    "lima_grp": "LImA_DATA",
                    "dataset_name": "data",
                    "compression": None,
                    "min_size":1,
                    "detector_name": "LImA Detector",
                    "metadata_grp": "detector_config"
                     }

    def __init__(self):
        """
        Constructor of an HDF5 writer:
        
        @param filename: name of the file
        @param hpath: name of the group: it will contain data (2-4D dataset), [tth|q|r] and pyFAI, group containing the configuration
        @param burst: exprected size of the dataset 
        """
        self._sem = threading.Semaphore()
        self._initialized = False
        self.filename = self.dirname = self.extension = self.subdir = self.hpath = self.lima_grp =None
        self.dataset_name = self.compression = self.min_size = self.detector_name = self.metadata_grp = None
        for kw, defval in self.CONFIG_ITEMS.items():
            self.__setattr__(kw, defval)

        self.hdf5 = None
        self.group = None
        self.dataset = None
        self.chunk = None
        self.shape = None

    def __repr__(self):
        out = ["HDF5 writer  %sinitialized" % ("" if self._initialized else "un")] + \
        ["%s: %s" % (k, self.__getattribute__(k)) for k in self.CONFIG_ITEMS]
        return os.linesep.join(out)

    def init(self, lima_cfg):
        """
        Initializes the HDF5 file for writing. Part of prepareAcq.
        
        @param lima_cfg: dictionnary with parameters coming from Lima at the "prepareAcq" 
        """
        with self._sem:

            self.filename = lima_cfg.get("directory")
            if not self.filename.endswith('/'): self.filename +=  '/'
            self.filename += lima_cfg.get("prefix")+'.h5'

            # silence non serious error messages, which are printed
            # because we use h5py in a new thread (not in the main one)
            # this is a bug seems to be fixed on newer version !!
            # see this h5py issue 206: https://code.google.com/p/h5py/issues/detail?id=206
            h5py._errors.silence_errors()

            try:
                self.hdf5 = h5py.File(self.filename, 'a')

            except IOError:
                os.unlink(self.filename)
                print ("here e e ")
                self.hdf5 = h5py.File(self.filename)


            prefix = lima_cfg.get("prefix") or self.CONFIG_ITEMS["hpath"]
            if not prefix.endswith("_"):
                prefix+="_"
            entries = len([i.startswith(prefix) for i in self.hdf5])
            self.hpath = posixpath.join("%s%04d"%(prefix,entries),self.lima_grp)
            self.group = self.hdf5.require_group(self.hpath)
            self.group.parent.attrs["NX_class"] = "NXentry"
            self.group.attrs["NX_class"] = "NXdata"
            cfg_grp = self.hdf5.require_group(posixpath.join(self.hpath, self.metadata_grp))
            cfg_grp["detector_name"] = numpy.string_(self.detector_name)
            for k, v in lima_cfg.items():
                if type(v) in types.StringTypes:
                    cfg_grp[k] = numpy.string_(v)
                else:
                    cfg_grp[k] = v
            self.number_of_frames = (max(1, lima_cfg["number_of_frames"]) + self.min_size - 1) // self.min_size * self.min_size
            self.min_size = max(1, self.min_size)
            self.shape = (self.number_of_frames , lima_cfg.get("dimY", 1), lima_cfg.get("dimX", 1))
            self.chunk = (self.min_size, lima_cfg.get("dimY", 1), lima_cfg.get("dimX", 1))
            if "dtype" in lima_cfg:
                self.dtype = numpy.dtype(lima_cfg["dtype"])
            else:
                self.dtype = numpy.int32
            self.dataset = self.group.require_dataset(self.dataset_name, self.shape, dtype=self.dtype, chunks=self.chunk,
                                                      maxshape=(None,) + self.chunk[1:])
            self.dataset.attrs["interpretation"] = "image"
            self.dataset.attrs["metadata"] = self.metadata_grp
            self.dataset.attrs["signal"] = "1"
            self.group.parent["title"] = numpy.string_("Raw frames")
            self.group.parent["program"] = numpy.string_("LImA HDF5 plugin")
            self.group.parent["start_time"] = numpy.string_(self.getIsoTime())

    def getIsoTime(self, forceTime=None):
        """
        @param forceTime: enforce a given time (current by default)
        @type forceTime: float
        @return: the current time as an ISO8601 string
        @rtype: string
        """
        if forceTime is None:
            forceTime = time.time()
        localtime = time.localtime(forceTime)
        gmtime = time.gmtime(forceTime)
        tz_h = localtime.tm_hour - gmtime.tm_hour
        tz_m = localtime.tm_min - gmtime.tm_min
        return "%s%+03i:%02i" % (time.strftime("%Y-%m-%dT%H:%M:%S", localtime), tz_h, tz_m)
    
    

    def flush(self):
        """
        Update some data like axis units and so on.        
        """
        with self._sem:
            if not self.hdf5:
                err = 'No opened file'
                raise RuntimeError(err)
            if "stop_time" in self.group.parent:
                del  self.group.parent["stop_time"]
            self.group.parent["stop_time"] = numpy.string_(self.getIsoTime())
            self.hdf5.flush()

    def close(self):
        self.flush()
        with self._sem:
            if self.hdf5:
                self.hdf5.close()
                self.hdf5 = None
                self.group = None
                self.dataset = None
                self.chunk = None
                self.size = None

    def write(self, data, index=0):
        """
        Minimalistic method to limit the overhead.
        """
        with self._sem:
            if index >= self.dataset.shape[0]:
                self.dataset.resize(index + 1, axis=0)
            self.dataset[index] = data

    def isLastFrame(self, index=0):
        return index == self.number_of_frames-1


class HDF5SavingDeviceServer(BasePostProcess) :
    HDF5_TASK_NAME = 'HDF5SavingTask'
    Core.DEB_CLASS(Core.DebModApplication, 'HDF5Saving')
    def __init__(self, cl, name):
        self.__HDF5SavingTask = None
        self.__HDF5Sink = None
        self.__HDF5Writer = None
        self.get_device_properties(self.get_device_class())
        self.__extension = ""
        self.__subdir = None
        BasePostProcess.__init__(self, cl, name)
        HDF5SavingDeviceServer.init_device(self)

    def set_state(self, state) :
        if(state == PyTango.DevState.OFF) :
            if(self.__HDF5SavingTask) :
                self.__HDF5Writer.close()
                self.__callback = None
                self.__HDF5Writer = None
                self.__HDF5SavingTask = None
                ctControl = _control_ref()
                extOpt = ctControl.externalOperation()
                extOpt.delOp(self.HDF5_TASK_NAME)
        elif(state == PyTango.DevState.ON) :
            if not self.__HDF5SavingTask:
                try:
                    ctControl = _control_ref()
                    extOpt = ctControl.externalOperation()
                    self.__HDF5SavingTask = extOpt.addOp(Core.USER_SINK_TASK,
                                                         self.HDF5_TASK_NAME,
                                                         self._runLevel)
                    if not self.__HDF5Sink:
                        self.__HDF5Writer = HDF5Writer()
                        self.__HDF5Sink = HDF5Sink(self.__HDF5Writer)
                    self.__HDF5SavingTask.setSinkTask(self.__HDF5Sink)
                    self.__callback = PrepareAcqCallback(ctControl, self.__HDF5Sink)
                    self.__HDF5SavingTask.registerCallback(self.__callback)

                except:
                    import traceback
                    traceback.print_exc()
                    return
        PyTango.Device_4Impl.set_state(self, state)


    def Reset(self) :
        if(self.__HDF5Sink) :
            self.__HDF5Sink.reset()

    def read_Parameters(self, the_att):
        """
        Called  when reading the "Parameters" attribute
        """
        if self.__HDF5Sink:
            the_att.set_value(self.__HDF5Sink.__repr__())
        else:
            the_att.set_value("No HDF5 Sink processlib active for the moment")


class HDF5SavingDeviceServerClass(PyTango.DeviceClass) :
        #        Class Properties
    class_property_list = {
        }


    #    Device Properties
    device_property_list = {
        }


    #    Command definitions
    cmd_list = {
        'Start':
        [[PyTango.DevVoid, ""],
         [PyTango.DevVoid, ""]],
        'Stop':
        [[PyTango.DevVoid, ""],
         [PyTango.DevVoid, ""]],
        }


    #    Attribute definitions
    attr_list = {
        'RunLevel':
            [[PyTango.DevLong,
            PyTango.SCALAR,
            PyTango.READ_WRITE]],
        'Parameters':
            [[PyTango.DevString,
            PyTango.SCALAR,
            PyTango.READ]],
        }
#------------------------------------------------------------------
#    HDF5SavingDeviceServerClass Constructor
#------------------------------------------------------------------
    def __init__(self, name):
        PyTango.DeviceClass.__init__(self, name)
        self.set_type(name)

_control_ref = None
def set_control_ref(control_class_ref):
    global _control_ref
    _control_ref = control_class_ref

def get_tango_specific_class_n_device() :
    return HDF5SavingDeviceServerClass, HDF5SavingDeviceServer

