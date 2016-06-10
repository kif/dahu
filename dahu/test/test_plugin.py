#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
from __future__ import with_statement, print_function

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "10/06/2016"
__status__ = "development"

import unittest
from . import utilstest
logger = utilstest.getLogger("test_plugin")
from dahu.plugin import Plugin, plugin_from_function
from dahu.factory import plugin_factory


class TestPlugin(unittest.TestCase):
    def test_plugin(self):
        "Test a stand alone (dummy-) plugin"
        p = Plugin()
        p.setup()
        p.process()
        p.teardown()
        logger.debug(p.output)

    def test_plugin_from_function(self):
        "Test a plugin from a function"
        def cube(a):
            return a * a * a
        plugin_name = plugin_from_function(cube)
        logger.debug(plugin_name)
        plugin = plugin_factory(plugin_name)
        # plugin from functions are callable:
        assert plugin(a=5) == 125


def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(TestPlugin("test_plugin"))
    testSuite.addTest(TestPlugin("test_plugin_from_function"))

    return testSuite

if __name__ == '__main__':
    mysuite = suite()
    runner = unittest.TextTestRunner()
    runner.run(mysuite)
