#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""Dahu: Data analysis server controlled from Tango 
"""
from __future__ import with_statement, print_function, absolute_import, division

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "19/10/2016"
__status__ = "production"


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
        return sloctime
    else:
        sloctime = time.strftime("%Y-%m-%dT%H:%M:%S", localtime)
        return "%s%+03i:%02i" % (sloctime, tz_h, tz_m)


def get_workdir(basedir=None):
    """
    Creates a working directory
    
    @param basedir: temporary directory
    @return: path of the working directory 
    """
    workdir = globals()["workdir"]
    basedir = basedir or ''
    if not(workdir) or not (workdir.startswith(basedir)):
        if not basedir:
            basedir = tempfile.gettempdir()
        foldername = os.path.basename(basedir)
        if foldername.startswith("dahu_") and len(foldername) == 24:
            # likely the time has already been added
            workdir = os.path.abspath(basedir)
        else:
            subdir = "dahu_%s" % get_isotime(for_path=True)
            workdir = os.path.join(basedir, subdir)
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
        globals()["workdir"] = workdir
    return workdir


def fully_qualified_name(obj):
    """
    Return the fully qualified name of an object
    
    @param obj: any python object
    @return: the full name as a string
    """
    if "__module__" not in dir(obj):  # or "__name__" not in dir(obj):
        obj = obj.__class__
    module = obj.__module__.lower()
    name = obj.__name__.lower()
    return module + "." + name
