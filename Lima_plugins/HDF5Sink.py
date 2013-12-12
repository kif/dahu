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
        lima_cfg = {"dimX":x,
                    "dimY":y,
                    "binX":binX,
                    "binY":binY}

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
        print("self._task._worker: %s" % self._task._worker)
        if (self._task._worker) is None :
            centerX = x // 2
            centerY = y // 2
            ai = pyFAI.AzimuthalIntegrator()
            ai.setFit2D(1000, centerX=centerX, centerY=centerY, pixelX=1, pixelY=1)
            worker = pyFAI.worker.Worker(ai)
            worker.unit = "r_mm"
            worker.method = "lut_ocl_gpu"
            worker.nbpt_azim = 360
            worker.nbpt_rad = 500
            worker.output = "numpy"
            print("Worker updated")
            self._task._worker = worker
        else:
            worker = self._task._worker
        worker.reconfig(shape=(y, x), sync=True)
        if self._task._writer:
            config = self._task._worker.get_config()
            self._task._writer.init(fai_cfg=config, lima_cfg=lima_cfg)
            self._task._writer.flush(worker.radial, worker.azimuthal)


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


    def setJsonConfig(self, jsonconfig):
        self._writer.setJsonConfig(jsonconfig)

class HDF5Writer(object):
    """
    Class allowing to write HDF5 Files.    
    """
    CONFIG_ITEMS = ["filename", "dirname", "extension", "subdir", "hpath"]
    DATASET_NAME = "data"
    def __init__(self, filename, hpath="data", burst=None):
        """
        Constructor of an HDF5 writer:
        
        @param filename: name of the file
        @param hpath: name of the group: it will contain data (2-4D dataset), [tth|q|r] and pyFAI, group containing the configuration
        @param burst: exprected size of the dataset 
        """
        self.filename = filename
        self._sem = threading.Semaphore()
        self.dirname = None
        self.subdir = None
        self.extension = extension
        self.fai_cfg = {}
        self.lima_cfg = {}
        self.hpath = hpath
        self.hdf5 = None
        self.group = None
        self.dataset = None
        self.chunk = None
        self.shape = None
        self.ndim = None

    def __repr__(self):
        return "HDF5 writer on file %s:%s %sinitialized" % (self.filename, self.hpath, "" if self._initialized else "un")

    def init(self, fai_cfg=None, lima_cfg=None):
        """
        Initializes the HDF5 file for writing
        @param fai_cfg: the configuration of the worker as a dictionary
        """
        Writer.init(self, fai_cfg, lima_cfg)
        with self._sem:
            #TODO: this is Debug statement
            open("fai_cfg.json", "w").write(json.dumps(self.fai_cfg, indent=4))
            open("lima_cfg.json", "w").write(json.dumps(self.lima_cfg, indent=4))
            self.fai_cfg["nbpt_rad"] = self.fai_cfg.get("nbpt_rad", 1000)
            if h5py:
                try:
                    self.hdf5 = h5py.File(self.filename)
                except IOError: #typically a corrupted HDF5 file !
                    os.unlink(self.filename)
                    self.hdf5 = h5py.File(self.filename)
            else:
                logger.error("No h5py library, no chance")
                raise RuntimeError("No h5py library, no chance")
            self.group = self.hdf5.require_group(self.hpath)
            self.group.attrs["NX_class"] = "NXentry"
            self.pyFAI_grp = self.hdf5.require_group(posixpath.join(self.hpath, self.CONFIG))
            self.pyFAI_grp.attrs["desc"] = "PyFAI worker configuration"
            for key, value in self.fai_cfg.items():
                if value is None:
                    continue
                try:
                    self.pyFAI_grp[key] = value
                except:
                    print("Unable to set %s: %s" % (key, value))
                    self.close()
                    sys.exit(1)
            rad_name, rad_unit = str(self.fai_cfg.get("unit", "2th_deg")).split("_", 1)
            self.radial_values = self.group.require_dataset(rad_name, (self.fai_cfg["nbpt_rad"],), numpy.float32)
            if self.fai_cfg.get("nbpt_azim", 0) > 1:
                self.azimuthal_values = self.group.require_dataset("chi", (self.fai_cfg["nbpt_azim"],), numpy.float32)
                self.azimuthal_values.attrs["unit"] = "deg"
                self.radial_values.attrs["interpretation"] = "scalar"
                self.radial_values.attrs["long name"] = "Azimuthal angle"

            self.radial_values.attrs["unit"] = rad_unit
            self.radial_values.attrs["interpretation"] = "scalar"
            self.radial_values.attrs["long name"] = "diffraction radial direction"
            if self.fast_scan_width:
                self.fast_motor = self.group.require_dataset("fast", (self.fast_scan_width,) , numpy.float32)
                self.fast_motor.attrs["long name"] = "Fast motor position"
                self.fast_motor.attrs["interpretation"] = "scalar"
                self.fast_motor.attrs["axis"] = "1"
                self.radial_values.attrs["axis"] = "2"
                if self.azimuthal_values is not None:
                    chunk = 1, self.fast_scan_width, self.fai_cfg["nbpt_azim"], self.fai_cfg["nbpt_rad"]
                    self.ndim = 4
                    self.azimuthal_values.attrs["axis"] = "3"
                else:
                    chunk = 1, self.fast_scan_width, self.fai_cfg["nbpt_rad"]
                    self.ndim = 3
            else:
                self.radial_values.attrs["axis"] = "1"
                if self.azimuthal_values is not None:
                    chunk = 1, self.fai_cfg["nbpt_azim"], self.fai_cfg["nbpt_rad"]
                    self.ndim = 3
                    self.azimuthal_values.attrs["axis"] = "2"
                else:
                    chunk = 1, self.fai_cfg["nbpt_rad"]
                    self.ndim = 2

            if self.DATASET_NAME in self.group:
                del self.group[self.DATASET_NAME]
            shape = list(chunk)
            if self.lima_cfg.get("number_of_frames", 0) > 0:
                if self.fast_scan_width is not None:
                    size[0] = 1 + self.lima_cfg["number_of_frames"] // self.fast_scan_width
                else:
                    size[0] = self.lima_cfg["number_of_frames"]
            self.dataset = self.group.require_dataset(self.DATASET_NAME, shape, dtype=numpy.float32, chunks=chunk,
                                                      maxshape=(None,) + chunk[1:])
            if self.fai_cfg.get("nbpt_azim", 0) > 1:
                self.dataset.attrs["interpretation"] = "image"
            else:
                self.dataset.attrs["interpretation"] = "spectrum"
            self.dataset.attrs["signal"] = "1"
            self.chunk = chunk
            self.shape = chunk
            name = "Mapping " if self.fast_scan_width else "Scanning "
            name += "2D" if self.fai_cfg.get("nbpt_azim", 0) > 1 else "1D"
            name += " experiment"
            self.group["title"] = name
            self.group["program"] = "PyFAI"
            self.group["start_time"] = getIsoTime()



    def flush(self, radial=None, azimuthal=None):
        """
        Update some data like axis units and so on.
        
        @param radial: position in radial direction
        @param  azimuthal: position in azimuthal direction
        """
        with self._sem:
            if not self.hdf5:
                raise RuntimeError('No opened file')
            if radial is not None:
                if radial.shape == self.radial_values.shape:
                    self.radial_values[:] = radial
                else:
                    logger.warning("Unable to assign radial axis position")
            if azimuthal is not None:
                if azimuthal.shape == self.azimuthal_values.shape:
                    self.azimuthal_values[:] = azimuthal
                else:
                    logger.warning("Unable to assign azimuthal axis position")
            self.hdf5.flush()

    def close(self):
        with self._sem:
            if self.hdf5:
                self.flush()
                self.hdf5.close()
                self.hdf5 = None

    def write(self, data, index=0):
        """
        Minimalistic method to limit the overhead.
        """
        with self._sem:
            if self.dataset is None:
                logger.warning("Writer not initialized !")
                return
            if self.azimuthal_values is None:
                data = data[:, 1] #take the second column only aka I
            if self.fast_scan_width:
                index0, index1 = (index // self.fast_scan_width, index % self.fast_scan_width)
                if index0 >= self.dataset.shape[0]:
                    self.dataset.resize(index0 + 1, axis=0)
                self.dataset[index0, index1] = data
            else:
                if index >= self.dataset.shape[0]:
                    self.dataset.resize(index + 1, axis=0)
                self.dataset[index] = data
