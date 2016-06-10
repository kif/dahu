#!/usr/bin/python
#coding: utf-8

# Stuff to read C216 time frame generator
from __future__ import with_statement, print_function, division

__authors__ = [ "Jérôme Kieffer"]
__contact__ = "jerome.kieffer@esrf.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20140425"
__status__ = "devel"

import os
if not os.environ.get("TANGO_HOST"):
    os.environ["TANGO_HOST"] = "saxs1new:20000"
import logging
logging.basicConfig()
logger = logging.getLogger("C216")
import PyTango
import numpy

class C216(object):
    """
    Instance to manage a C216 device server
    """
    def __init__(self, device_id):
        self.id = device_id
        self.device = PyTango.DeviceProxy(device_id)
        self.last_status = {}

    def is_running(self):
        """
        Returns True if the device is running 
        """
        res = self.device.CompStatus("Tango::RUNNING")
        logger.debug("device %s in state %s" % (self.id, res))
        return res == "Tango::ON"

    @property
    def frames(self):
        """returns the number of available frames"""
        return self.status["FRAMES"]

    @property
    def status(self):
        """
        read the status of the device
        """
        status = self.device.GetCompleteConfigAndStatus()
        res = {"TFU_MODE":          status[0],
               "TFU_FRAME_MODE":    status[1],
               "START_MODE":        status[2],
               "FRAMES"     :       status[14] // 2,
               "CYCLES"      :      status[17],
               "CHANNELS"    :      status[19],
               "ACQ_BANK"       :   status[65],
               "TFU_CONFIG":        status[:28],
               "DIGI_IN":           status[28:44],
               "ANA_IN":            status[44:60]}
        self.last_status = res
        return res

    def read_scalers(self):
        """
        Return a 2D array of 
        """
        if self.is_running():
            logger.warning("C216 needs to be stopped before reading scalers")
#            self.stop
        nb_frames = int(self.frames)
        data = self.device.ReadScalersForNLiveFrames([0, nb_frames - 1])
        data.shape = self.frames, -1
        return data

if __name__ == "__main__":
    c216 = C216("id02/c216/0")

    print(c216.is_running())
    print(c216.read_scalers())
    print(c216.last_status)

