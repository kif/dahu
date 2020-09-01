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

import unittest
from . import utilstest
logger = utilstest.getLogger("test_plugin")
from ..plugin import Plugin, plugin_from_function
from ..factory import plugin_factory
from ..job import Job


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

    def test_wait_for(self):
        "Test the synchronization of plugins"
        p = Plugin()
        print(dir(p))
        print(p.__class__.__module__)
        p.wait_for(42) #this job does not exist, fails with a warning:
        
        # Test synchonization with finished job
        j = Job("example.square", {"x": 5})
        j.start()
        p.wait_for(j.id)
        
        #Test failure when it does not start (timeout)
        p.TIMEOUT=0.5
        j = Job("example.square", {"x": 6})
        try:
            p.wait_for(j.id)
        except RuntimeError as err:
            logger.debug("Failed as expected with: %s", err)
        else:
            raise RuntimeError("Expected to fail !")
        
        
        

        


def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(TestPlugin("test_plugin"))
    testSuite.addTest(TestPlugin("test_plugin_from_function"))
    testSuite.addTest(TestPlugin("test_wait_for"))

    return testSuite

if __name__ == '__main__':
    mysuite = suite()
    runner = unittest.TextTestRunner()
    runner.run(mysuite)
