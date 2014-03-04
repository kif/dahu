#!/usr/bin/env python
# -*- coding: utf8 -*-
#
from __future__ import with_statement, print_function

__doc__ = """Data Analysis Highly tailored for Upbl09a 
            """
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20140303"
__status__ = "development"
version = "0.1"

import time
import os
import tempfile
workdir = None

def get_isotime(forceTime=None, for_path=False):
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
    if for_path:
        sloctime = time.strftime("%Y-%m-%dT%Hh%Mm%S", localtime)
        return "%s%+03i%02i" % (sloctime, tz_h, tz_m)
    else:
        sloctime = time.strftime("%Y-%m-%dT%H:%M:%S", localtime)
        return "%s%+03i:%02i" % (sloctime, tz_h, tz_m)

def get_workdir(basedir=""):
    """
    Creates a working directory
    """
    workdir = globals()["workdir"]
    if not(workdir) or not (workdir.startswith(basedir)):
        if not basedir:
            basedir = tempfile.gettempdir()
        subdir = "dahu_%s" % get_isotime(for_path=True)
        workdir = os.path.join(basedir, subdir)
        os.makedirs(workdir)
        globals()["workdir"] = workdir
    return workdir
