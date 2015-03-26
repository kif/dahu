#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Data Analysis plugin for focus analysis 
"""

from __future__ import with_statement, print_function
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "26/03/2015"
__status__ = "development"
version = "0.0.0"

import os
import numpy
from dahu.plugin import Plugin
from dahu.factory import register
from threading import Semaphore
import pyFAI
import logging
logger = logging.getLogger("plugin.focus")
import scipy, scipy.misc, scipy.stats


@register
class Plugin(Plugin):
    """Plugin in charge of calculating a descriptor of the focus quality related to the decay of the power spectrum of an image
    

    Structure of the input data:
input = {
        "filename":"toto.jpg"
        }

"""

    _dictGeo = {}  # key:shape, value= ai
    _sem = Semaphore()

    def process(self):
        img = scipy.misc.imread(self.input["filename"], True)
        shape = img.shape
        nr = numpy.sqrt(shape[0] * shape[0] + shape[1] * shape[1]) // 2
        if shape not in self._dictGeo:
            with self._sem:
                if shape not in self._dictGeo:
                    ai = pyFAI.AzimuthalIntegrator()
                    ai.setFit2D(1000, shape[0] // 2, shape[1] // 2, pixelX=10, pixelY=10)
                    self._dictGeo[shape] = ai
        ai = self._dictGeo[shape]
        ft = numpy.fft.fft2(img)
        fts = numpy.fft.fftshift(ft)
        r, i = ai.integrate1d(abs(fts) ** 2, nr, unit="r_mm")
        value = -scipy.stats.linregress(numpy.log(r), numpy.log(i))[0]
        self.output["filename"] = self.input["filename"]
        self.output["value"] = value
