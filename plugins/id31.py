#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import with_statement, print_function

"""Data Analysis plugin tailored for ID31

* integrate_simple: simple demo of a  
* integrate: a more advances 
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "24/03/2016"
__status__ = "development"
version = "0.1.0"

import os
import numpy
from dahu.plugin import Plugin, plugin_from_function
from dahu.factory import register
from threading import Semaphore
import logging
logger = logging.getLogger("plugin.pyFAI")

try:
    import pyFAI
except ImportError:
    logger.error("Failed to import PyFAI: download and install it from pypi")
try:
    import fabio
except ImportError:
    logger.error("Failed to import Fabio: download and install it from pypi")


def integrate_simple(poni_file, image_file, curve_file, nbins=1000):
    """Simple azimuthal integration for a single frame (very inefficient)
    
    :param poni_file: configuration of the geometry
    :param image_file: 
    :param curve_file: output file
    :param nbins: number of output bins
    """
    ai = pyFAI.load(poni_file)
    img = fabio.open(image_file).data
    ai.integrate1d(img, nbins, filename=curve_file, unit="2th_deg", method="splitpixel")
    return {"out_file": curve_file}
plugin_from_function(integrate_simple)


# Use thre register decorator to make it available from Dahu
@register
class PluginIntegrate(Plugin):
    """This is the basic plugin of PyFAI for azimuthal integration
    
    Input parameters:
    :param poni_file: configuration of the geometry
    :param input_files:
    :param 
    """
    _dictGeo = {}  # key:tuple(ai.param), value= ai
    _sem = Semaphore()

    def __init__(self):
        """
        """
        Plugin.__init__(self)
        self.shared = None
        self.strOutputFile = None
        self.ai = None  # this is the azimuthal integrator to use
        self.data = None
        self.mask = None
        self.dummy = None
        self.delta_dummy = None
        self.dest_dir = None

    def setup(self, kwargs):
        Plugin.setup(self, kwargs)
        logger.debug("PluginIntegrate.setup")
        Plugin.setup(self, kwargs)
        if "output_dir" not in self.input:
            self.log_error("output_dir not in input")
        self.dest_dir = os.path.abspath(self.input["output_dir"])
        ai = pyFAI.AzimuthalIntegrator()


        ########################################################################
        # Choose the azimuthal integrator
        ########################################################################

        with self.__class__._sem:
            if tuple(ai.param) in self.__class__._dictGeo:
                self.ai = self.__class__._dictGeo[tuple(ai.param)]
            else:
                self.__class__._dictGeo[tuple(ai.param)] = ai
                self.ai = ai


    def process(self):
        Plugin.process(self)
        logger.debug("PluginIntegrate.process")
        data = EDUtilsArray.getArray(self.dataInput.input)
        if self.dataInput.saxsWaxs and self.dataInput.saxsWaxs.value.lower().startswith("s"):
            out = self.ai.saxs(self.data,
                               nbPt=self.nbPt,
                               filename=self.strOutputFile,
                               mask=self.mask,
                               dummy=self.dummy,
                               delta_dummy=self.delta_dummy)
        else:
            out = self.ai.xrpd(self.data,
                               nbPt=self.nbPt,
                               filename=self.strOutputFile,
                               mask=self.mask,
                               dummy=self.dummy,
                               delta_dummy=self.delta_dummy)
        self.npaOut = numpy.hstack((i.reshape(-1, 1) for i in out if i is not None))

    def teardown(self):
        Plugin.teardown(self)
        logger.debug("PluginPyFAIv1_0.teardown")
        # Create some output data
