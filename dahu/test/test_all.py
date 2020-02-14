#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""Test suite for all dahu modules."""

from __future__ import with_statement, print_function

__authors__ = ["Jérôme Kieffer"]
__contact__ = "jerome.kieffer@esrf.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__data__ = "10/06/2016"

import sys
import unittest
from .utilstest import getLogger
logger = getLogger(__file__)

from . import test_job
from . import test_plugin


def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_job.suite())
    testSuite.addTest(test_plugin.suite())
    return testSuite

if __name__ == '__main__':
    mysuite = suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
