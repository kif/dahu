"""
X-ray photon correlation spectroscopy plugin for ID02 

"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "25/03/2020"
__status__ = "development"
__version__ = "0.1.0"

import logging
logger = logging.getLogger("id02.xpcs")

from dahu.plugin import Plugin
import numpy
import h5py
try:
    import numexpr
except ImportError as err:
    logger.error("NumExpr module is missing, some calculation will be slower")
else:
    numexpr = None
import fabio
from pyFAI.detectors import Detector
from pyFAI.geometry import Geometry


class xpcs(Plugin):
    """This plugin does pixel correlation for XPCS and averages the signal from various bins provided in the qmask. 

Minimalistic example:
{
    "data_file": "Janus_Eiger500k_raw.h5",
    "dark_file": None,
    "result_file": "Janus_Eiger500k_xpcs.h5",
    "sample": {
        "name": "FAB_PPG",
        "composition": "FAB_PPG425_250",
        "temperature": "300K"
    },
    "detector": {
        "name": "eiger500k",
        "pixel": 75e-6,
        "mask": None,
        "flatfield": None
    },
    "experiment_setup":{
        "geometry": "SAXS",
        "wavelength": 1.538, #Å ?
        "detector_distance": 2.3, #m?
        "lagtime": 0.1, #?
        "firstq": 0.0036,
        "widthq": 4e-4,
        "stepq": 1e-4,
        "numberq": 27,
        "q_mask": "qmask.npy",
        "beamstop_mask": "mask.npy" ,
        "directbeam_x": 104, #pixel
        "directbeam_y": 157, #pixel 
    },
    "correlator":{
        "name": "MatMulCorrelator",
        "method": "intensity"
        #"low_threshold": 50
        #"bottom_ADU": 1420
        #"top_ADU": 2200
        #"photon_ADU"": 1970
        #"max_number"": 150000
        #"ttcf": 12
    
    }
             
}
"""

    def __init__(self):
        Plugin.__init__(self)
        self.dataset
        self.shape = None  # shape of every image
        self.nframes = None  # Number of input frames
        self.qmask = None  # contains the numpy array with the qmask

    def process(self, kwargs=None):
        Plugin.process(self, kwargs)
        #

    def make_qmask(self):
        "create the q_mask from the geometry"

        experiment_setup = self.input.get("experiment_setup", {})
        if experiment_setup.get("q_mask"):
            qmask = fabio.open(experiment_setup["q_mask"]).data
        else:
            detector_section = self.input.get("detector", {})
            pixel_size = detector_section.get("pixel")
            if pixel_size is None:
                self.log_error("Pixel size is mandatory in detector description section")
            detector = Detector(pixel1=pixel_size, pixel2=pixel_size, max_shape=self.shape)
            wavelength = experiment_setup.get("wavelength")
            if wavelength is None:
                self.log_error("wavelength is mandatory in experiment_setup section")
            else:
                wavelength *= 1e-10  # Convert Å in m
            distance = experiment_setup.get("detector_distance")
            if distance is None:
                self.log_error("detector_distance is mandatory in experiment_setup section")
            directbeam_x = experiment_setup.get("directbeam_x")
            directbeam_y = experiment_setup.get("directbeam_y")
            if (directbeam_x is None) or (directbeam_y is None):
                self.log_error("directbeam_[xy] is mandatory in experiment_setup section")

            geometry = Geometry(distance, directbeam_y * pixel_size, directbeam_x * pixel_size,
                                detector=detector, wavelength=wavelength)
            q_array = geometry.center_array(self.shape, unit="q_nm^-1")

            firstq = experiment_setup.get("firstq", 0)
            widthq = experiment_setup.get("widthq")
            stepq = experiment_setup.get("stepq", 0)
            numberq = experiment_setup.get("numberq", (1 << 16) - 1)  # we plan to store the qmask as uint16
            if widthq is None:
                self.log_error("widthq is mandatory in experiment_setup section")
            if numexpr is not None:
                
            else:
                qmaskf = (q_array - firstq) / (widthq + stepq)
                qmaskf[qmaskf<0] = 0
                qmaskf[qmaskf>(numberq+1)] = 0
                qmaskf[(qmaskf%1)>widthq/(widthq + stepq)] = 0
                qmask_numpy = qmaskf.astype(dtype=numpy.uint16)

                self.log_warning("numexpr is missing, calculation slower")

        return qmask

    def teardown(self):
        self.output["result_file"] = self.result_filename
        if self.nxs:
            self.nxs.close()
        Plugin.teardown(self)
