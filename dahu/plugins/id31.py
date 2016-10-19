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
__date__ = "19/10/2016"
__status__ = "development"
version = "0.1.0"

import os
import numpy
from ..plugin import Plugin, plugin_from_function
from ..factory import register
from ..cache import DataCache
from threading import Semaphore
import logging
logger = logging.getLogger("plugin.pyFAI")
import json

try:
    import pyFAI
    from pyFAI.worker import make_ai
except ImportError:
    logger.error("Failed to import PyFAI: download and install it from pypi")
try:
    import fabio
except ImportError:
    logger.error("Failed to import Fabio: download and install it from pypi")
try:
    import calculate_flux as flux
except ImportError:
    logger.error("Failed to import flux calculator for ID31")


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



# Use the register decorator to make it available from Dahu

@register
class Integrate(Plugin):
    """This is the basic plugin of PyFAI for azimuthal integration
    
    Input parameters:
    :param poni_file: configuration of the geometry
    :param input_files:
    :param 
    
    Typical JSON file:
    {"poni_file": "/tmp/example.poni",
     "input_files": ["/tmp/file1.edf", "/tmp/file2.edf"],
     "monitor_values": [1, 1.1],
     "npt": 2000,
     "unit": "2th_deg",
    }
    """
    _ais = DataCache()  # key: str(a), value= ai

    def __init__(self):
        """
        """
        Plugin.__init__(self)
        self.ai = None  # this is the azimuthal integrator to use
        self.dest_dir = None
        self.json_data = None
        self.ntp = 3000
        self.input_files = []
        self.method = "full_ocl_csr"
        self.unit = "q_nm^-1"
        self.output_files = []
        self.mask = ""
        self.wavelength = None
        self.dummy = -1
        self.delta_dummy = 0
        self.polarization_factor = None
        self.do_SA = False
        self.norm = 1e12

    def setup(self, kwargs):
        logger.debug("Integrate.setup")
        Plugin.setup(self, kwargs)

        if "output_dir" not in self.input:
            self.log_error("output_dir not in input")
        # this needs to be added in the SPEC macro
        self.dest_dir = os.path.abspath(self.input["output_dir"])
        if "json" not in self.input:
            self.log_error("json not in input")
        json_path = self.input.get("json", "")
        if not os.path.exists(json_path):
            self.log_error("Integration setup file (JSON): %s does not exist" % json_path, do_raise=True)
        self.json_data = json.load(open(json_path))

        ai = make_ai(self.json_data)
        stored = self._ais.get(str(ai), ai)
        if stored is ai:
            self.ai = stored
        else:
            self.ai = stored.__deepcopy__()

        self.npt = int(self.json_data.get("npt", self.npt))
        self.unit = self.json_data.get("unit", self.unit)
        self.wavelength = self.json_data.get("wavelength", self.wavelength)
        if os.path.exists(self.json_data["mask"]):
            self.mask = self.json_data.get("mask", self.mask)
        self.dummy = self.json_data.get("val_dummy", self.dummy)
        self.delta_dummy = self.json_data.get("delta_dummy", self.delta_dummy)
        if self.json_data["do_polarziation"]:
            self.polarization_factor = self.json_data.get("polarization_factor", self.polarization_factor)
        self.do_SA = self.json_data.get("do_SA", self.do_SA)
        self.norm = self.json_data.get("norm", self.norm)  # need to be added in the spec macro

    def process(self):
        Plugin.process(self)
        logger.debug("Integrate.process")
        for fname in self.input_files:
            if not os.path.exists(fname):
                self.log_error("image file: %s does not exist, skipping" % fname,
                               do_raise=False)
                continue
            basename = os.path.splitext(os.path.basename(fname))[0]
            destination = os.path.join(self.dest_dir, basename + ".dat")
            fimg = fabio.open(fname)
            if self.wavelength is not None:
                monitor = self.getMon(fimg.header, self.wavelength) / self.norm
            else:
                monitor = 1.0
            self.ai.integrate1d(fimg.data, npt=self.npt, method=self.method,
                                safe=False,
                                filename=destination,
                                normalization_factor=monitor,
                                unit=self.unit,
                                dummy=self.dummy,
                                delta_dummy=self.delta_dummy,
                                polarization_factor=self.polarization_factor,
                                correctSolidAngle=self.do_SA
                                )
            self.output_files.append(destination)

    def teardown(self):
        Plugin.teardown(self)
        logger.debug("Integrate.teardown")
        # Create some output data
        self.output["output_files"] = self.output_files

    @staticmethod
    def getMon(header, lam):
        strCount = header['counter_mne'].split()
        strCountPos = header['counter_pos'].split()
        E = 4.13566766225e-15 * 299792458 / lam / 1000
        return flux.main(float(strCountPos[strCount.index('mondio')]), E)

