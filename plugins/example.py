#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example of plugins made from functions 
"""
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "2021-05-10"
__status__ = "development"
version = "0.2.0"

import time
from dahu.plugin import plugin_from_function, Plugin
from dahu.factory import register
import logging
logger = logging.getLogger("plugin.example")


@register
class Cube(Plugin):

    def process(self):
        Plugin.process(self)
        if self.input is None:
            logger.warning("input is None")
        x = self.input.get("x", 0)
        self.output["result"] = x * x * x


def square(x):
    return x * x


plugin_from_function(square)


def noop(x):
    return x


plugin_from_function(noop)


@register
class Sleep(Plugin):
    """
    Plugin hat sleeps for a while ...
    """

    def process(self):
        if self.input is None:
            logger.warning("input is None, defaulting to 10 sec")
            t = 10
        else:
            t = self.input.get("t", 10)
        time.sleep(t)
