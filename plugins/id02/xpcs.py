"""
X-ray photon correlation spectroscopy plugin for ID02 

"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "24/06/2020"
__status__ = "development"
__version__ = "0.1.0"

import os
import json
import logging
logger = logging.getLogger("id02.xpcs")

from dahu.plugin import Plugin
import numpy
try:
    import numexpr
except ImportError as err:
    logger.error("NumExpr module is missing, some calculation will be slower")
    numexpr = None

import hdf5plugin
import h5py
import fabio
from .common import Nexus, get_isotime
from pyFAI.detectors import Detector
from pyFAI.geometry import Geometry
from pyFAI.units import CONST_hc
from dynamix import version as dynamix_version
from dynamix.correlator import dense
# Dummy factory for correlators
CORRELATORS = {i: getattr(dense, i)
               for i in dir(dense)
               if i.endswith("Correlator")}
COMPRESSION = hdf5plugin.Bitshuffle()


class XPCS(Plugin):
    """This plugin does pixel correlation for XPCS and averages the signal from various bins provided in the qmask. 

Minimalistic example:
{
    "plugin_name": "id02.xpcs",
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
        "mask": null,
        "flatfield": null
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
        self.result_filename = None
        self.timestep = None

    def process(self):
        Plugin.process(self)
        self.dataset = self.read_data()
        self.nframes = self.dataset.shape[0]
        self.shape = self.dataset.shape[1:]
        self.qmask = self.make_qmask()
        Correlator = self.get_correlator()
        correlator = Correlator(self.shape, self.nframes, qmask=self.qmask)
        results = correlator.correlate(self.dataset[...])
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
        if "measurement" in entry:
            measurement = entry["measurement"]
        else:
                self.log_error("No measurement in entry: %s of data_file: %s" % (entry, data_file))
        if "data" in measurement:
            dataset = measurement["data"]
        else:
            self.log_error("No dataset in measurement: %s of data_file: %s" % (entry, data_file))
            dataset = None

        instruments = nxs.get_class(entry, "NXinstrument")
        if len(instruments) == 0:
            self.log_error("No NXinstrument in entry: %s of data_file: %s" % (entry, data_file))
        instrument = instruments[0]
        detectors = nxs.get_class(instrument, "NXdetector")
        if len(detectors) != 1:
            self.log_error("No NXdetector in entry: %s of data_file: %s" % (entry, data_file))
        detector = detectors[0]
        header = detector["header"]
        zero = numpy.array([0])
        self.timestep = float(header.get("acq_expo_time", zero)[()]) + float(header.get("acq_latency_time", zero)[()])
        return dataset

    def make_qmask(self):
        "create the q_mask from the geometry"

        experiment_setup = self.input.get("experiment_setup", {})
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
        self.ai = geometry

        firstq = experiment_setup.get("firstq", 0)
        widthq = experiment_setup.get("widthq", 0)
        stepq = experiment_setup.get("stepq", 0)
        numberq = experiment_setup.get("numberq", (1 << 16) - 2)  # we plan to store the qmask as uint16

        if experiment_setup.get("q_mask"):
            qmask = fabio.open(experiment_setup["q_mask"]).data
        else:
            q_array = geometry.center_array(self.shape, unit=self.unit)

            detector_maskfile = detector_section.get("mask", '')
            if os.path.exists(detector_maskfile):
                detector_mask = fabio.open(detector_maskfile).data
            else:
                detector_mask = detector.mask
                if detector_mask is None:
                    detector_mask = numpy.zeros(self.shape, dtype=numpy.int8)
            beamstop_maskfile = experiment_setup.get("beamstop_mask", "")
            if os.path.exists(beamstop_maskfile, ""):
                beamstop_mask = fabio.open(beamstop_maskfile).data
            else:
                beamstop_mask = numpy.zeros(self.shape, dtype=numpy.int8)
            mask = numpy.logical_or(detector_mask, beamstop_mask)

            if widthq is None:
                self.log_error("widthq is mandatory in experiment_setup section")
            if numexpr is None:
                qmaskf = (q_array - firstq) / (widthq + stepq)
                qmaskf[qmaskf < 0] = 0
                qmaskf[qmaskf > (numberq + 1)] = 0
                qmaskf[(qmaskf % 1) > widthq / (widthq + stepq)] = 0
                qmaskf[mask] = 0
                qmask = qmaskf.astype(dtype=numpy.uint16)
                self.log_warning("numexpr is missing, calculation is slower")
            else:
                qmaskf = numexpr.evaluate("(q_array - firstq) / (widthq + stepq)")
                qmask = numexpr.evaluate("where(qmaskf<0, 0, where(qmaskf>(numberq+1),0, where((qmaskf%1)>(widthq/(widthq + stepq)), 0, where(mask, 0, qmaskf))))",
                         out=numpy.empty(q_array.shape, dtype=numpy.uint16),
                         casting="unsafe")

        self.qrange = firstq + widthq / 2.0 + numpy.arange(qmask.max()) * (widthq + stepq)

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

        if os.path.exists(result_file):
            os.unlink(result_file)

        with Nexus(result_file, mode="w", creator="dahu") as nxs:
            entry_grp = nxs.new_entry(entry="entry",
                                      program_name=self.input.get("plugin_name", "dahu"),
                                      title="XPCS experiment",
                                      force_time=self.start_time)
            nxs.h5.attrs["default"] = entry_grp.name

            # Sample description, provided by the input
            sample_grp = nxs.new_class(entry_grp, "sample", "NXsample")
            for key, value in self.input.get("sample", {}).items():
                sample_grp[key] = value

            # Process 0: measurement group
            measurement_grp = nxs.new_class(entry_grp, "0_measurement", "NXdata")
            measurement_grp["images"] = h5py.ExternalLink(os.path.relpath(os.path.abspath(self.dataset.file.filename),
                                                                          os.path.dirname(os.path.abspath(result_file))),
                                                          self.dataset.name)
            measurement_grp.attrs["signal"] = "images"

            # Instrument
            instrument_grp = nxs.new_instrument(entry_grp, "ID02")
            instrument_grp["name"] = "TruSAXS"
            source_grp = nxs.new_class(instrument_grp, "ESRF", "NXsource")
            source_grp["type"] = "Synchrotron X-ray source"
            source_grp["name"] = "European Synchrotron Radiation Facility"
            source_grp["probe"] = "X-ray"
            monochromator_grp = nxs.new_class(instrument_grp, "monochromator", "NXmonochromator")
            monochromator_grp["description"] = "Cryogenically cooled Si-111 channel-cut monochromator"
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
            crystal_grp = nxs.new_class(monochromator_grp, "Si-111", "NXcrystal")
            crystal_grp["usage"] = "Bragg"
            crystal_grp["type"] = "Si"
            crystal_grp["temperature"] = numpy.float32(77.0)
            crystal_grp["temperature"].attrs["unit"] = "K"
            crystal_grp["reflection"] = numpy.ones(3, dtype=numpy.int8)
            # Storage of the result
            xpcs_grp = nxs.new_class(entry_grp, "1_XPCS", "NXprocess")
            xpcs_grp["sequence_index"] = 1
            xpcs_grp["program"] = "dynamix"
            xpcs_grp["version"] = dynamix_version
            xpcs_grp["date"] = get_isotime()
            xpcs_grp.create_dataset("correlator", data=self.correlator_name)

            config_grp = nxs.new_class(xpcs_grp, "configuration", "NXnote")
            config_grp["type"] = "text/json"
            config_grp["data"] = json.dumps(self.input, indent=2, separators=(",\r\n", ": "))

            xpcs_data = nxs.new_class(xpcs_grp, "results", "NXdata")
            qmask_ds = xpcs_data.create_dataset("qmask", self.qmask.shape, chunks=self.qmask.shape, **COMPRESSION)
            qmask_ds[...] = self.qmask
            qmask_ds.attrs["interpretation"] = "image"
            qmask_ds.attrs["long_name"] = "mask with bins averaged (0=masked-out)"

            entry_grp.attrs["default"] = xpcs_grp.attrs["default"] = xpcs_data.name
            result_ds = xpcs_data.create_dataset("g2", result.shape, chunks=result.shape, **COMPRESSION)
            result_ds[...] = result
            result_ds.attrs["interpretation"] = "spectrum"
            qrange_ds = xpcs_data.create_dataset("q", data=self.qrange)
            qrange_ds.attrs["interpretation"] = "scalar"
            qrange_ds.attrs["unit"] = self.unit
            qrange_ds.attrs["long_name"] = "Scattering vector q (%s)" % self.unit

            trange = numpy.arange(result.shape[-1]) * self.timestep
            trange_ds = xpcs_data.create_dataset("t", data=trange)
            trange_ds.attrs["interpretation"] = "scalar"
            trange_ds.attrs["unit"] = "s"
            trange_ds.attrs["long_name"] = "Time t (s)"

            xpcs_data.attrs["signal"] = "g2"
            xpcs_data.attrs["axes"] = ["q", "t"]
            xpcs_data["title"] = "g₂(q, t)"
        self.result_filename = result_file

    def teardown(self):
        if self.result_filename:
            self.output["result_file"] = self.result_filename
        try:
            self.dataset.file.close()
        except Exception as err:
            self.log_warning("%s Unable to close dataset file: %s" % (type(err), err))
        self.qmask = None
        self.ai = None
        self.correlator_name = None
        Plugin.teardown(self)
