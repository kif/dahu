"""
Common code used by all plugins from ID02
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "03/09/2020"
__status__ = "development"
__version__ = "1.0.0"

import sys
import logging
import h5py
import pyFAI.io
pyFAI.io.logger.setLevel(logging.ERROR)
if pyFAI.version_info < (0, 20):
    from .nexus import Nexus, get_isotime
else:
    from pyFAI.io import Nexus, get_isotime

# silence non serious error messages, which are printed
# because we use h5py in a new thread (not in the main one)
# this is a bug seems to be fixed on newer version !!
# https://github.com/h5py/h5py/issues/206

try:
    h5py._errors.silence_errors()
except:
    pass

StringTypes = (str, bytes)

def ensure_str(junk):
    "return a unicode string, regardless to the input"
    if isinstance(junk, bytes):
        return junk.decode()
    else:
        return str(junk)
    
""" Remaining code from former single file plugin 
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

Plugins for ID02:

* Distortion correction
* Metadata saving (C216)
* single detector processing


from __future__ import with_statement, print_function, division

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "16/03/2020"
__status__ = "production"
version = "0.8"

import h5py
import logging
import numpy
import os
import posixpath
import pyFAI
try:
    import bitshuffle
    import bitshuffle.h5
except:
    bitshuffle = None
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
import pyFAI.worker
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

from dahu.job import Job
from dahu.cache import DataCache


if sys.version_info[0] < 3:
    StringTypes = (str, unicode)
else:
    StringTypes = (str, bytes)

logger = logging.getLogger("dahu.id02")
# set loglevel at least at INFO
if logger.getEffectiveLevel() > logging.INFO:
    logger.setLevel(logging.INFO)
if logger.getEffectiveLevel() < logging.root.level:
    logger.setLevel(logging.root.level)
    
"""
