#!/usr/bin/env python3
#

"""Test suite for all dahu modules."""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "jerome.kieffer@esrf.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__data__ = "10/06/2016"

import sys
import unittest

from . import test_cache, test_job, test_plugin
from .utilstest import getLogger

logger = getLogger(__file__)


def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_job.suite())
    testSuite.addTest(test_plugin.suite())
    testSuite.addTest(test_cache.suite())
    return testSuite

if __name__ == '__main__':
    mysuite = suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
