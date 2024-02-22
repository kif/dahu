#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "22/02/2022"
__status__ = "production"

import unittest
from . import utilstest
logger = utilstest.getLogger(__name__)
from ..cache import DataCache


class TestCache(unittest.TestCase):
    def test_cache(self):       
        b0 = DataCache(5, borg=True)
        b1 = DataCache(7, borg=True)
        self.assertEqual(b0.max_size, 5, "borg initialized")
        self.assertEqual(b1.max_size, 5, "borg not re-initialized")
        
        n0 = DataCache(5, borg=False)
        n1 = DataCache(7, borg=False)
        self.assertEqual(n0.max_size, 5, "Normal class behaviour")
        self.assertEqual(n1.max_size, 7, "Normal class behaviour")

        b0["Iam"] = "borg"
        n0["Iam"] = "normal"
        self.assertEqual(b0.get("Iam"), "borg", "straight forwards")
        self.assertEqual(b1.get("Iam"), "borg", "borg behavour")
        self.assertEqual(n0.get("Iam"), "normal", "Normal class`")
        self.assertEqual(n1.get("Iam"), None, "Normal class`")
        
        for i in range(6):
            b0[i] = i
            n0[i] = i
        
        self.assertEqual(b0.get("Iam"), None, "object dropped")
        self.assertEqual(b1.get("Iam"), None, "object dropped, borg")
        self.assertEqual(n0.get("Iam"), None, "object dropped")
        self.assertEqual(n1.get("Iam"), None, "Normal class`")
        
def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(TestCache("test_cache"))
    return testSuite

if __name__ == '__main__':
    mysuite = suite()
    runner = unittest.TextTestRunner()
    runner.run(mysuite)
