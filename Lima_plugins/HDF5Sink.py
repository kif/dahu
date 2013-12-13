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

import os, json, distutils.util, sys, threading, logging
logger = logging.getLogger("lima.hdf5")
# set loglevel at least at INFO
if logger.getEffectiveLevel() > logging.INFO:
    logger.setLevel(logging.INFO)
import numpy
from Lima import Core
#from Utils import BasePostProcess
import sys
import os, json, distutils.util
from os.path import dirname
from pyFAI.io import getIsoTime
import h5py

class StartAcqCallback(Core.SoftCallback):
    """
    Class managing the connection from a 
    Lima.Core.CtControl.prepareAcq() to the configuration of the various tasks
    
    Example of usage:
    cam = Basler.Camera(ip)
    iface = Basler.Interface(cam)
    ctrl = Core.CtControl(iface)
    processLink = LinkPyFAI(worker, writer)
    extMgr = ctrl.externalOperation()
    myOp = self.extMgr.addOp(Core.USER_LINK_TASK, "pyFAILink", 0)
    myOp.setLinkTask(processLink)
    callback = StartAcqCallback(ctrl, processLink)
    myOp.registerCallback(callback)
    acq.setAcqNbFrames(0)
    acq.setAcqExpoTime(1.0)
    ctrl.prepareAcq() #Configuration called here !!!!
    ctrl.startAcq()

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
                    "rotation":im.getRotation(),
                    "mode": im.getMode(),
                    "image_type": im.getImageType()}
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
        lima_cfg["number_of_frames"] = acq.getAcqNbFrames() #to check.
        lima_cfg["exposure_time"] = acq.getAcqExpoTime()
        #ROI see: https://github.com/esrf-bliss/Lima/blob/master/control/include/CtAcquisition.h
        if self._task._writer:
            self._task._writer.init(lima_cfg=lima_cfg)
            self._task._writer.flush()


class HDF5Sink(Core.Processlib.SinkTaskBase):
    """
    This is a ProcessLib task which is a sink: 
    it saves the image into a HDF5 stack.
    
    """

    def __init__(self, worker=None, writer=None):
        Core.Processlib.SinkTaskBase.__init__(self)
        self._config = {}
        self._writer = writer
        if writer  is None:
            logger.error("Without a writer, SinkPyFAI will just dump all data")

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
        rData = Core.Processlib.Data()
        rData.frameNumber = data.frameNumber
        rData.buffer = self._worker.process(data.buffer)
        if self._writer: #optional HDF5 writer
            self._writer.write(rData.buffer, rData.frameNumber)


    def setConfig(self, config_dict=None):
        """
        Set the "static" configuration like filename and so on.
        
        @param config_dict: dict or json-serialized dict or file containing this dict.
        """
        self._writer.setConfig(jsonconfig)

class HDF5Writer(object):
    """
    Class allowing to write HDF5 Files.    
    """
    CONFIG_ITEMS = {"filename": None,
                    "dirname":None,
                    "extension": ".h5",
                    "subdir":None,
                    "hpath": "raw",
                    "dataset_name": "data",
                    "compression": None,
                    "min_size":1,
                    "detector_name": "LImA Detector",
                     }

    def __init__(self, **config):
        """
        Constructor of an HDF5 writer:
        
        @param filename: name of the file
        @param hpath: name of the group: it will contain data (2-4D dataset), [tth|q|r] and pyFAI, group containing the configuration
        @param burst: exprected size of the dataset 
        """
        self._sem = threading.Semaphore()
        self._initialized = False
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

    def init(self, lima_cfg=None):
        """
        Initializes the HDF5 file for writing. Part of prepareAcq.
        
        @param lima_cfg: dictionnary with parameters coming from Lima at the "prepareAcq" 
        """
        with self._sem:
            if h5py:
                try:
                    self.hdf5 = h5py.File(self.filename)
                except IOError:
                    logger.error("typically a corrupted HDF5 file ! : %s" % self.filename)
                    os.unlink(self.filename)
                    self.hdf5 = h5py.File(self.filename)
            else:
                err = "No h5py library, no chance"
                logger.error(err)
                raise RuntimeError(err)
            self.group = self.hdf5.require_group(self.hpath)
            self.group.attrs["NX_class"] = "NXentry"
            cfg_grp = self.hdf5.require_group(posixpath.join(self.hpath, "detector_config"))
            cfg_grp["detector_name"] = self.detector_name
            for k, v in lima_cfg.items:
                cfg_grp[k] = v
            #TODO: calculate chunks, size and datatype ...
            self.dataset = self.group.require_dataset(self.dataset_name, self.shape, dtype=numpy.float32, chunks=chunk,
                                                      maxshape=(None,) + chunk[1:])
            self.dataset.attrs["interpretation"] = "image"
            self.dataset.attrs["signal"] = "1"
            self.chunk = chunk
            self.shape = chunk
            self.group["title"] = "Raw frames"
            self.group["program"] = "LImA"
            self.group["start_time"] = getIsoTime()



    def flush(self):
        """
        Update some data like axis units and so on.
        
        @param radial: position in radial direction
        @param  azimuthal: position in azimuthal direction
        """
        with self._sem:
            if not self.hdf5:
                err = 'No opened file'
                logger.error(err)
                raise RuntimeError(err)
            self.group["stop_time"] = getIsoTime()
            self.hdf5.flush()

    def close(self):
        with self._sem:
            if self.hdf5:
                self.flush()
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
            if self.dataset is None:
                logger.warning("Writer not initialized !")
                return
            if index >= self.dataset.shape[0]:
                self.dataset.resize(index + 1, axis=0)
            self.dataset[index] = data

    def setConfig(self, config_dict=None):
        """
        Set the "static" configuration like filename and so on.
        
        @param config_dict: dict or JSON-serialized dict or file containing this dict.
        """

        if type(config_dict) in types.StringTypes:
            if os.path.isfile(config_dict):
                config = json.load(open(config_dict, "r"))
            else:
                 config = json.loads(config_dict)
        else:
            config = dict(config_dict)
        for k, v in  config.items():
            if k in self.CONFIG_ITEMS:
                self.__setattr__(k, v)
