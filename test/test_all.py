#!/usr/bin/env python
# -*- coding: utf8 -*-
#
from __future__ import with_statement, print_function
__doc__ = """
Test suite for all pyFAI modules.
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "jerome.kieffer@esrf.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__data__ = "20140324"

import sys
import unittest
from utilstest import UtilsTest, getLogger
logger = getLogger(__file__)

from test_job import test_suite_all_job
from test_plugin import test_suite_all_plugin

def test_suite_all():
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_suite_all_job())
    testSuite.addTest(test_suite_all_plugin())
    return testSuite

if __name__ == '__main__':
    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

