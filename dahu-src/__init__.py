#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Data Analysis Highly tailored for Upbl09a 
"""

from __future__ import with_statement, print_function
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "18/03/2015"
__status__ = "development"

from ._version import version, version_info, hexversion, date

import sys, logging
logging.basicConfig()

from . import utils
from . import factory
from . import plugin
from . import job
