#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
from __future__ import with_statement, print_function

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "01/09/2020"
__status__ = "development"

import os
import unittest
from . import utilstest
logger = utilstest.getLogger(__name__)
from .. import job


class TestJob(unittest.TestCase):
    def test_plugin_from_function(self):
        self.called = False
        dico = {"plugin_name": "example.square",
                "x": 5}
        j = job.Job("example.square", dico)
        j.connect_callback(self.callback)
        logger.info(j)
        j.start()
        j.join()
        logger.info(j)
        if "error" in j.output_data:
            logger.error(os.linesep.join(j.output_data["error"]))

        logger.info(j.input_data)
        self.assert_(self.called)
        print(j.output_data)
        self.assertEqual(j.output_data["result"], 25, "result OK")

    def callback(self, *args, **kwargs):
        logger.info("callback actually called with  arguments %s and kwargs %s" % (args, kwargs))
        assert len(args) == 1
        self.called = True


def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(TestJob("test_plugin_from_function"))
    return testSuite

if __name__ == '__main__':
    mysuite = suite()
    runner = unittest.TextTestRunner()
    runner.run(mysuite)
