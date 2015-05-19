#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Data Analysis Highly tailored for Upbl09a 
"""

from __future__ import with_statement, print_function, absolute_import, division

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "19/05/2015"
__status__ = "production"

from ._version import version, version_info, hexversion, date

import sys, logging
logging.basicConfig()

from . import utils
from . import factory
from . import plugin
from . import job
