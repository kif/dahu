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

import os
import json
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
import hdf5plugin
import fabio
from pyFAI.detectors import Detector
from pyFAI.geometry import Geometry
from pyFAI.io import Nexus, get_isotime
from pyFAI.units import CONST_hc
from dynamix import version as dynamix_version
from dynamix.correlator import dense
# Dummy factory for correlators
CORRELATORS = {(i, getattr(dense, i))
               for i in dir(dense)
               if i.endswith("Correlator")}
COMPRESSION = hdf5plugin.Bitshuffle()


class xpcs(Plugin):
    """This plugin does pixel correlation for XPCS and averages the signal from various bins provided in the qmask. 

Minimalistic example:
{
    "data_file": "Janus_Eiger500k_raw.h5",
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
        "wavelength": 1.538, #Å ?
        "detector_distance": 2.3, #m?
        "unit": "q_nm^-1",
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
    }
             
}
"""

    def __init__(self):
        Plugin.__init__(self)
        self.start_time = get_isotime()
        self.nxs = None  # output data file
        self.dataset = None
        self.shape = None  # shape of every image
        self.nframes = None  # Number of input frames
        self.qmask = None  # contains the numpy array with the qmask
        self.unit = "q_nm^-1"
        self.ai = None
        self.correlator_name = None

    def process(self, kwargs=None):
        Plugin.process(self, kwargs)
        self.dataset = self.read_data()
        self.nframes = self.dataset.shape[0]
        self.shape = self.dataset.shape[1:]
        self.qmask = self.make_qmask()
        Correlator = self.get_correlator()
        correlator = Correlator(self.shape, self.nframes, qmask=self.qmask)
        results = correlator.correlate(self.frames[...])
        self.save_results(results)

    def read_data(self):
        data_file = self.input.get("data_file")
        if data_file is None:
            self.log_error("data_file not provided, mandatory")
        elif not os.path.exists(data_file):
            self.log_error("data_file does not exist")
        nxs = Nexus(data_file, mode="r")
        entries = nxs.get_entries()
        if len(entries) == 0:
            self.log_error("No entry in data_file %s" % data_file)
        entry = entries[0]
        measurments = nxs.get_class(entry, "NXmeasurment")
        if len(measurments) == 0:
            self.log_error("No NXmeasurment in entry: %s of data_file: %s" % (entry, data_file))
        measurment = measurments[0]
        if "data" in measurment:
            dataset = measurment["data"]
        else:
            self.log_error("No dataset in measurment: %s of data_file: %s" % (entry, data_file))
            dataset = None
        return dataset

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

            self.unit = experiment_setup.get("unit", self.unit)

            geometry = Geometry(distance, directbeam_y * pixel_size, directbeam_x * pixel_size,
                                detector=detector, wavelength=wavelength)
            q_array = geometry.center_array(self.shape, unit=self.unit)

            firstq = experiment_setup.get("firstq", 0)
            widthq = experiment_setup.get("widthq")
            stepq = experiment_setup.get("stepq", 0)
            numberq = experiment_setup.get("numberq", (1 << 16) - 2)  # we plan to store the qmask as uint16

            # TODO: manage the different masks!

            if widthq is None:
                self.log_error("widthq is mandatory in experiment_setup section")
            if numexpr is not None:
                qmaskf = numexpr.evaluate("(q_array - firstq) / (widthq + stepq)")
                qmask = numexpr.evaluate("where(qmaskf<0, 0, where(qmaskf>(numberq+1),0, where((qmaskf%1)>(widthq/(widthq + stepq)), 0, qmaskf)))",
                         out=numpy.empty(q_array.shape, dtype=numpy.uint16),
                         casting="unsafe")
            else:
                qmaskf = (q_array - firstq) / (widthq + stepq)
                qmaskf[qmaskf < 0] = 0
                qmaskf[qmaskf > (numberq + 1)] = 0
                qmaskf[(qmaskf % 1) > widthq / (widthq + stepq)] = 0
                qmask = qmaskf.astype(dtype=numpy.uint16)
                self.log_warning("numexpr is missing, calculation is slower")
        self.ai = geometry
        return qmask

    def get_correlator(self):
        correlator_section = self.input.get("correlator", {})
        correlator_name = correlator_section.get("name")
        if correlator_name is None:
            self.log_error("No correlator name in input")
        if correlator_name not in CORRELATORS:
            self.log_error("Correlator requested %s is not part of the available ones: %s" % (correlator_name, ", ".join(CORRELATORS.keys())))
        self.correlator_name = correlator_name
        return CORRELATORS[correlator_name]

    def save_results(self, result):
        "The structure is inspired from the pyFAI data format which stores reduced data"
        result_file = self.input.get("result_file")
        if result_file is None:
            input_file = os.path.abspath(self.dataset.file.filename)
            a, b = os.path.splitext(input_file)
            result_file = a + "_xpcs" + b
            self.log_warning("No destination file provided, saving in %s" % result_file)
        with Nexus(result_file, mode="w") as nxs:
            entry_grp = nxs.new_entry(self, entry="entry",
                                      program_name=self.input.get("plugin_name", "dahu"),
                                      title="XPCS experiment",
                                      force_time=self.start_time)
            nxs.h5.attrs["default"] = entry_grp.name

            # Sample description, provided by the input
            sample_grp = nxs.new_class(entry_grp, self.sample.name, "NXsample")
            for key, value in self.input.get("sample", {}):
                sample_grp[key] = value

            # Process 0: Measurement group
            measurement_grp = nxs.new_class(entry_grp, "0_measurement", "NXdata")
            measurement_grp["images"] = h5py.ExternalLink(self.dataset.file.filename, self.dataset.name)
            measurement_grp.attrs["signal"] = "images"

            # Instrument
            instrument_grp = nxs.new_instrument(entry_grp, "ID02")
            instrument_grp["name"] = "TruSAXS"
            source_grp = nxs.new_class(instrument_grp, "ESRF", "NXsource")
            source_grp["radiation"] = "Synchrotron X-ray source"
            source_grp["name"] = "European Synchrotron Radiation Facility"
            source_grp["probe"] = "X-ray"
#             current = numpy.ascontiguousarray(self.input.get("storage_ring_current", []), dtype=numpy.float32)
#             current_ds = source_grp.create_dataset("current",
#                                                    data=current)
#             current_ds.attrs["units"] = "mA"
#             current_ds.attrs["interpretation"] = "spectrum"
            monochromator_grp = nxs.new_class(instrument_grp, "DCM", "NXmonochromator")
            monochromator_grp["description"] = "Cryogenically cooled Si-111 double-crystal monochromator "
            wl = self.ai.wavelength * 1e10
            resolution = 2e-4
            wl_ds = monochromator_grp.create_dataset("wavelength", data=numpy.float32(wl))
            wl_ds.attrs["units"] = "Å"
            wl_ds.attrs["resolution"] = resolution
            wle_ds = monochromator_grp.create_dataset("wavelength_error", data=numpy.float32(wl * resolution))
            wle_ds.attrs["units"] = "Å"
            nrj_ds = monochromator_grp.create_dataset("energy", data=CONST_hc / wl)
            nrj_ds.attrs["units"] = "keV"
            nrj_ds.attrs["resolution"] = resolution
            nrje_ds = monochromator_grp.create_dataset("energy_error", data=CONST_hc / wl * resolution)
            nrje_ds.attrs["units"] = "keV"

            # Storage of the result
            xpcs_grp = nxs.new_class(entry_grp, "1_XPCS", "NXprocess")
            xpcs_grp["sequence_index"] = 1
            xpcs_grp["program"] = "dynamix"
            xpcs_grp["version"] = dynamix_version
            xpcs_grp["date"] = get_isotime()
            xpcs_grp.create_dataset("correlator", self.correlator_name)
            xpcs_grp.create_dataset("correlator_config", data=json.dumps({}))  # TODO
            xpcs_grp.create_dataset("qmask", data=self.qmask, **COMPRESSION)
            xpcs_grp.attrs["interpretation"] = "image"
            xpcs_data = nxs.new_class(xpcs_grp, "results", "NXdata")
            xpcs_grp.attrs["default"] = xpcs_data.name
            result_ds = xpcs_data.create_dataset("data", data=result, **COMPRESSION)
            result_ds.attrs["interpretation"] = "spectrum"

    def teardown(self):
        self.output["result_file"] = self.result_filename
        try:
            self.dataset.file.close()
        except Exception as err:
            self.log_warning("%s Unable to close dataset file: %s" % (type(err), err))
        self.qmask = None
        self.ai = None
        self.correlator_name = None
        Plugin.teardown(self)
