#!/usr/bin/env python
# -*- coding: utf8 -*-
#
from __future__ import with_statement, print_function

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20140324"
__status__ = "development"

import sys, os, unittest
import utilstest
logger = utilstest.getLogger("test_job")
import dahu.job

class TestJob(unittest.TestCase):
    def test_plugin_from_function(self):
        dico = {"plugin_name": "example.square", "x": 5 }
        j = dahu.job.Job(dico)
        print(j)
        j.start()
        print(j.output_data)


def test_suite_all_job():
    testSuite = unittest.TestSuite()
    testSuite.addTest(TestJob("test_plugin_from_function"))

    return testSuite

if __name__ == '__main__':
    mysuite = test_suite_all_job()
    runner = unittest.TextTestRunner()
    runner.run(mysuite)
