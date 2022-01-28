"""
Distortion correction and azimuthal integration for a single detector 

"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "28/01/2022"
__status__ = "development"
__version__ = "0.9.2"

import os
import threading
import shutil
import posixpath
from collections import namedtuple
import json
import logging
logger = logging.getLogger("id02.single_detector")
import numpy
import numexpr
import h5py
from dahu import version as dahu_version
from dahu.plugin import Plugin
from dahu.utils import get_isotime, fully_qualified_name
from dahu.cache import DataCache
import fabio
import pyFAI, pyFAI.utils
from pyFAI.worker import Worker, DistortionWorker, PixelwiseWorker
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from .common import StringTypes, Nexus, ensure_str
try:
    import hdf5plugin
except ImportError:
    COMPRESSION = {}
else:
    COMPRESSION = hdf5plugin.Bitshuffle()

numexpr.set_num_threads(8)

CacheKey = namedtuple("CacheKey", ["ai", "shape"])


class CacheValues:
    __slots__ = 'array', 'engine'

    def __init__(self, array):
        self.array = array
        self.engine = {}

    def update(self, array=None, engine=None):
        if array:
            self.array.update(array)
        if engine:
            self.engine.update(engine)


class SingleDetector(Plugin):
    """This plugin does all processing needed for a single camera

Minimalistic example:
              {
 "npt1_rad": 1000,
 "c216_filename": "/nobackup/lid02gpu11/metadata/test.h5",
 "npt2_rad": 500,
 "DetectorName": "rayonix",
 "npt2_azim": 360,
 "to_save": "raw sub ave azim dist flat ave_log",
 "output_dir": "/nobackup/lid02gpu12/output",
 "WaveLength": 9.95058e-11,
 "image_file": "/nobackup/lid02gpu11/FRELON/test_laurent_saxs_0000.h5",
 "dark_filename": "/nobackup/lid02gpu11/FRELON/test_laurent_dark_0000.h5",
}

All possible options are :
Mandatory options:
"image_file": HDF5 file containing the raw data
"output_dir": directory where final data are saved
"to_save": list of processing to be performed like: "raw sub flat dist norm azim ave"
"c216_filename": File containing the metadata (I1, timimgs ...)
"npt2_azim": number on bins in q direction for azim file
"npt2_rad": number on bins in azimuthal direction for azim file
"npt1_rad": number of bins in q direction for ave file

Optional parameters:
"DetectorName": "rayonix",
"PSize_1": 2.4e-05,
"PSize_2": 2.4e-05,
"BSize_1":1,
"BSize_2":1,
"Center_1":512,
"Center_2":512,
"DDummy":1,
"SampleDistance":14.9522,
"WaveLength": 9.95058e-11,
"Dummy":-10,
"metadata_job": int: wait for this job to finish and retrieve metadata from it.
"regrouping_mask_filename": file containing the mask to be used to integration
"dark_filename" HDF5 file with dark data in it
'dark_filter"  can be 'quantil' to average between given lower and upper quantils
"dark_filter_quantil" ask for the median id 0.5, min if 0 or max if 1
"dark_filter_quantil_lower" lower quantil for averaging
"dark_filter_quantil_upper" upper quantil for averaging
"distortion_filename" Spline file with the distortion correction
"flat_filename" flat field for intensity response correction
"metadata_job": number of the metadata job to wait for (unused)
"scaling_factor": float (default=1) by which all intensity will be multiplied
"correct_solid_angle": True by default, set to 0/False to disable such correction
"correct_I1": True by default, set to false to deactivate scaling by Exposure time / transmitted intensity
"unit": "q_nm^-1" can be changed to "log(q)_m" for log(q) units
"variance_formula": calculate the variance from a formula like '0.1*(data-dark)+1.0' "

Unused and automatically generated:
"plugin_name':"id02.singledetector',
"job_id': 29,

Defined but yet unused keywords:
"RasterOrientation": int,
"Rot": float

Possible values for to_save:
----------------------------

* sub: subtract dark signal
* flat: subtract dark then divide by flat
* solid: subtract dark then divide by flat and absolute solid angle
* dist: subtract dark then divide by flat and absolute solid angle, finally rebin on a regular grid
* norm: subtract dark then divide by flat and absolute solid angle, finally rebin on a regular grid and divide intensity by "I1"
* azim: subtract dark then divide by flat and absolute solid angle, finally rebin on a regular q/Chi grid and divide intensity by "I1"
* ave: subtract dark then divide by flat and absolute solid angle, finally rebin on a regular q grid and divide intensity by "I1"

    """
    KEY_CONV = {"BSize": int,
                "Offset": int,
                "RasterOrientation": int,
                "Center": float,
                "DDummy": float,
                "Dummy": float,
                "PSize": float,
                "SampleDistance": float,
                "Wavelength": float,  # both are possible ?
                "WaveLength": float,
                "Rot": float
                }

    KEYS = ("BSize_1", "BSize_2", "Center_1", "Center_2", "DDummy", "DetectorName",
            "Dummy", "Offset_1", "Offset_2", "PSize_1", "PSize_2",
            "Rot_1", "Rot_2", "Rot_3",
            "RasterOrientation", "SampleDistance", "SaxsDataVersion", "Title", "WaveLength")
    TIMEOUT = 10
    cache = DataCache(5)
    REPROCESS_IGNORE = ["metadata_job"]

    def __init__(self):
        Plugin.__init__(self)
        self.ai = None
        self.distortion = None
        self.workers = {}
        self.output_ds = {}  # output datasets
        self.dest = None  # output directory
        self.I1 = None  # beam stop diode values
        self.t = None  # time of the start of the frame. Same shape as self.I1
        self.nframes = None
        self.to_save = ["raw", "ave"]  # by default only raw image and averaged one is saved
        self.input_nxs = None
        self.metadata_nxs = None
        self.images_ds = None
        self.metadata_plugin = None
        self.metadata = {}
        self.npt1_rad = None
        self.npt2_rad = None
        self.npt2_azim = None
        self.dark = None
        self.dark_filename = None
        self.flat_filename = None
        self.flat = None
        self.mask_filename = None
        self.distortion_mask = None
        self.regrouping_mask = None
        self.distortion_filename = None
        self.output_hdf5 = {}
        self.dist = 1.0
        self.absolute_solid_angle = None
        self.in_shape = None
        self.scaling_factor = 1.0
        self.correct_solid_angle = True
        self.correct_I1 = True
        self.dummy = None
        self.delta_dummy = None
        self.unit = "q_nm^-1"
        self.polarization = None
        self.cache_ai = None
        self.cache_dis = None
        self.variance_formula = None
        self.variance_function = lambda data, dark: None

    def setup(self, kwargs=None):
        """
        see class documentation
        """
        Plugin.setup(self, kwargs)
        if "output_dir" not in self.input:
            self.log_error("output_dir not in input")
        self.dest = os.path.abspath(self.input["output_dir"])

        if "unit" in self.input:
            self.unit = self.input.get("unit")

        if "metadata_job" in self.input:
            job_id = int(self.input.get("metadata_job", 0))
            self.wait_for(job_id)

        if not os.path.isdir(self.dest):
            os.makedirs(self.dest)
        c216_filename = os.path.abspath(self.input.get("c216_filename", ""))

        if (os.path.dirname(c216_filename) != self.dest) and (os.path.basename(c216_filename) not in os.listdir(self.dest)):
            self.output_hdf5["metadata"] = os.path.join(self.dest, os.path.basename(c216_filename))
            m = threading.Thread(target=shutil.copy, name="copy metadata", args=(c216_filename, self.dest))
            m.start()

        if "to_save" in self.input:
            to_save = self.input["to_save"][:]
            if type(to_save) in StringTypes:
                # fix a bug from spec ...
                self.to_save = [i.strip('[\\] ",') for i in to_save.split()]
                self.log_warning("processing planned: " + " ".join(self.to_save))
            else:
                self.to_save = to_save
        if "image_file" not in self.input:
            self.log_error("image_file not in input")
        self.image_file = self.input["image_file"]
        if not os.path.exists(self.image_file):
            if not self.image_file.startswith("/"):
                # prepend the dirname of the c216
                image_file = os.path.join(os.path.dirname(c216_filename), self.image_file)
                if os.path.exists(image_file):
                    self.image_file = image_file
                else:
                    self.log_error("image_file %s does not exist" % self.image_file)

        self.dark_filename = self.input.get("dark_filename")
        if "raw" in self.to_save:
            if os.path.dirname(self.image_file) != self.dest:
                t = threading.Thread(target=shutil.copy, name="copy raw", args=(self.image_file, self.dest))
                t.start()
            self.output_hdf5["raw"] = os.path.join(self.dest, os.path.basename(self.image_file))
            if type(self.dark_filename) in StringTypes and os.path.exists(self.dark_filename):
                if os.path.dirname(self.dark_filename) != self.dest:
                    d = threading.Thread(target=shutil.copy, name="copy dark", args=(self.dark_filename, self.dest))
                    d.start()
                self.output_hdf5["dark"] = os.path.join(self.dest, os.path.basename(self.dark_filename))
        self.scaling_factor = float(self.input.get("scaling_factor", 1.0))
        self.correct_solid_angle = bool(self.input.get("correct_solid_angle", True))
        self.correct_I1 = bool(self.input.get("correct_I1", True))
        self.I1, self.t = self.load_I1_t(c216_filename)

        # Variance formula: calculation of the function
        if "variance_formula" in self.input:
            self.variance_formula = self.input.get("variance_formula")
            if self.variance_formula:
                self.variance_function = numexpr.NumExpr(self.variance_formula,
                                                         [("data", numpy.float64),
                                                          ("dark", numpy.float64)])

    def process(self):
        dummy = None
        self.metadata = self.parse_image_file()
        shape = self.in_shape[-2:]

        if self.I1 is None:
            self.I1 = numpy.ones(shape, dtype=float)
        elif self.I1.size < self.in_shape[0]:
            ones = numpy.ones(self.in_shape[0], dtype=float)
            ones[:self.I1.size] = self.I1
            self.I1 = ones
        # update all metadata with the one provided by input
        for key, value in self.input.items():
            if key in self.KEYS:
                self.metadata[key] = value
        forai = {}
        for key in ("BSize_1", "BSize_2", "Center_1", "Center_2",
                    "PSize_1", "PSize_2", "Rot_1", "Rot_2", "Rot_3",
                    "SampleDistance", "WaveLength"):
            forai[key] = self.metadata.get(key)

        self.dist = self.metadata.get("SampleDistance")

        # Some safety checks, use-input are sometimes False !
        if self.dist < 0:
            # which is impossible
            self.log_warning(f"Found negative distance: {self.dist}, considering its absolute value")
            self.dist = forai["SampleDistance"] = abs(self.metadata.get("SampleDistance"))

        if forai["WaveLength"] < 0:
            self.log_warning(f"Found negative wavelength: {forai['WaveLength']}, considering its absolute value")
            forai["WaveLength"] = abs(self.metadata.get("WaveLength"))
        self.dummy = self.metadata.get("Dummy", self.dummy)
        self.delta_dummy = self.metadata.get("DDummy", self.delta_dummy)

        self.ai = AzimuthalIntegrator()
        self.ai.setSPD(**forai)

        self.distortion_filename = self.input.get("distortion_filename") or None
        if type(self.distortion_filename) in StringTypes:
            detector = pyFAI.detector_factory("Frelon", {"splineFile": self.distortion_filename})
            if tuple(detector.shape) != shape:
                self.log_warning("Binning needed for spline ? detector claims %s but image size is %s" % (detector.shape, shape))
                detector.guess_binning(shape)
            self.ai.detector = detector
        else:
            self.ai.detector.shape = self.in_shape[-2:]

        self.log_warning("AI:%s" % self.ai)

        self.cache_ai = CacheKey(str(self.ai), shape)

        if self.cache_ai in self.cache:
            cached_ai = self.cache.get(self.cache_ai)
            self.ai._cached_array.update(cached_ai.array)
            self.ai.engines.update(cached_ai.engine)
        else:
            # Initialize geometry:
            self.ai.qArray(shape)
            self.ai.chiArray(shape)
            self.ai.deltaQ(shape)
            self.ai.deltaChi(shape)
            self.ai.solidAngleArray(shape)
            # Store cached array in cache
            self.cache.get(self.cache_ai, CacheValues({})).update(array=self.ai._cached_array)

        if self.input.get("do_polarization"):
            self.polarization = self.ai.polarization(factor=self.input.get("polarization_factor"),
                                                     axis_offset=self.input.get("polarization_axis_offset", 0))

        # Static mask is used for distortion correction
        static_mask = self.ai.detector.mask
        if static_mask is None:
            static_mask = numpy.zeros(shape, dtype=bool)

        # Read and Process dark
        if isinstance(self.dark_filename, StringTypes) and os.path.exists(self.dark_filename):
            dark = self.read_data(self.dark_filename)
            if dark.ndim == 3:
                method = self.input.get("dark_filter")
                if method.startswith("quantil"):
                    lower = self.input.get("dark_filter_quantil_lower", 0)
                    upper = self.input.get("dark_filter_quantil_upper", 1)
                    self.dark = pyFAI.average.average_dark(dark, center_method=method, quantiles=(lower, upper))
                else:
                    if method is None:
                        method = "median"
                    self.dark = pyFAI.average.average_dark(dark, center_method=method)
            else:
                self.dark = dark
        elif type(self.dark_filename) in (int, float):
            self.dark = float(self.dark_filename)
        if numpy.isscalar(self.dark):
            self.dark = numpy.ones(self.ai.detector.shape) * self.dark
        self.ai.detector.set_darkcurrent(self.dark)

        # Read and Process Flat
        self.flat_filename = self.input.get("flat_filename")
        if type(self.flat_filename) in StringTypes and os.path.exists(self.flat_filename):
            if self.flat_filename.endswith(".h5") or self.flat_filename.endswith(".nxs") or self.flat_filename.endswith(".hdf5"):
                flat = self.read_data(self.flat_filename)
            else:
                flat_fabio = fabio.open(self.flat_filename)
                flat = flat_fabio.data
                dummy = flat_fabio.header.get("Dummy")
                try:
                    dummy = float(dummy)
                except:
                    self.log_error("Dummy value in mask is unconsistent %s" % dummy)
                    dummy = None
                ddummy = flat_fabio.header.get("DDummy")
                try:
                    ddummy = float(ddummy)
                except:
                    self.log_error("DDummy value in mask is unconsitent %s" % ddummy)
                    ddummy = 0

            if flat.ndim == 3:
                self.flat = pyFAI.average.average_dark(flat, center_method="median")
            else:
                self.flat = flat
            if (self.flat is not None) and (self.flat.shape != shape):
                binning = [j / i for i, j in zip(shape, self.flat.shape)]
                if tuple(binning) != (1, 1):
                    self.log_warning("Binning for flat is %s" % binning)
                    if max(binning) > 1:
                        binning = [int(i) for i in binning]
                        self.flat = pyFAI.utils.binning(self.flat, binsize=binning, norm=False)
                    else:
                        binning = [i // j for i, j in zip(shape, self.flat.shape)]
                        self.flat = pyFAI.utils.unBinning(self.flat, binsize=binning, norm=False)
            if numpy.isscalar(self.flat):
                self.flat = numpy.ones(shape) * self.flat
            self.ai.detector.set_flatfield(self.flat)
            # extend the static mask with dummy pixels from the flat-field image
            if dummy:
                if ddummy:
                    numpy.logical_or(static_mask, abs(self.flat - dummy) <= ddummy, out=static_mask)
                else:
                    numpy.logical_or(static_mask, self.flat == dummy, out=static_mask)
            self.distortion_mask = static_mask
        self.ai.detector.mask = self.distortion_mask

        # Read and Process mask for integration
        self.mask_filename = self.input.get("regrouping_mask_filename")
        if isinstance(self.mask_filename, StringTypes) and os.path.exists(self.mask_filename):
            try:
                mask_fabio = fabio.open(self.mask_filename)
            except:
                local_mask = self.read_data(self.mask_filename) != 0
            else:  # this is very ID02 specific !!!!
                dummy = mask_fabio.header.get("Dummy")
                try:
                    dummy = float(dummy)
                except:
                    self.log_error("Dummy value in mask is unconsitent %s" % dummy)
                    dummy = None
                ddummy = mask_fabio.header.get("DDummy")
                try:
                    ddummy = float(ddummy)
                except:
                    self.log_error("DDummy value in mask is unconsitent %s" % ddummy)
                    ddummy = 0
                if ddummy:
                    local_mask = abs(mask_fabio.data - dummy) <= ddummy
                else:
                    local_mask = mask_fabio.data == dummy
                self.dummy = dummy
                self.delta_dummy = ddummy
            if local_mask.ndim == 3:
                local_mask = pyFAI.average.average_dark(local_mask, center_method="median")
            if (local_mask is not None) and (local_mask.shape != shape):
                binning = [j / i for i, j in zip(shape, local_mask.shape)]
                if tuple(binning) != (1, 1):
                    self.log_warning("Binning for mask is %s" % binning)
                    if max(binning) > 1:
                        binning = [int(i) for i in binning]
                        local_mask = pyFAI.utils.binning(local_mask, binsize=binning, norm=True) > 0
                    else:
                        binning = [i // j for i, j in zip(shape, local_mask.shape)]
                        local_mask = pyFAI.utils.unBinning(local_mask, binsize=binning, norm=False) > 0
            self.regrouping_mask = numpy.logical_or(self.distortion_mask, local_mask)
            # self.log_warning("found %s pixel masked out" % (local_mask.sum()))
        self.ai.detector.mask = self.regrouping_mask
        # bug specific to ID02, dummy=0 means no dummy !
        if self.dummy == 0:
            self.dummy = None

        self.create_hdf5()
        self.process_images()

    def load_I1_t(self, mfile, correct_shutter_closing_time=True):
        """
        load the I1 data and timstamp for frame start from a metadata HDF5 file


        /entry_0001/id02/MCS/I1

        TODO: handle correction or not for shutter opening/closing time

        @param mfile: metadata HDF5 file
        @param correct_shutter_closing_time: set to true for integrating detector (CCD) and false for counting detector (Pilatus)
        @return: 2-tuple of array with I1 and t
        """
        if ("I1" in self.input):
            return numpy.array(self.input["I1"]), None

        if not os.path.exists(mfile):
            self.log_error("Metadata file %s does not exist" % mfile, do_raise=True)

        self.metadata_nxs = Nexus(mfile, "r")
        I1 = t = None
        if correct_shutter_closing_time:
            key = "Intensity1ShutCor"
        else:
            key = "Intensity1UnCor"
        for entry in self.metadata_nxs.get_entries():
            for instrument in self.metadata_nxs.get_class(entry, "NXinstrument"):
                if "MCS" in instrument:
                    mcs = instrument["MCS"]
                    if key in mcs:
                        I1 = numpy.array(mcs[key])
                if "TFG" in instrument:
                    tfg = instrument["TFG"]
                    if "delta_time" in tfg:
                        t = numpy.array(tfg["delta_time"])
                if (t is None) or (I1 is None):
                    I1 = t = None
                else:
                    break
            if (t is not None) and (I1 is not None):
                return I1, t
        return I1, t

    def parse_image_file(self):
        self.input_nxs = Nexus(self.image_file, "r")
        creator = self.input_nxs.h5.attrs.get("creator")
        if creator is None:
            return self.parse_image_file_old()
        elif creator.startswith("LIMA-"):
            version = creator.split("-")[1]
            # test on version of lima
            if version < "1.":
                self.log_warning("Suspicious version of LIMA: %s" % creator)
            return self.parse_image_file_lima()

    def parse_image_file_old(self):
        """Historical version, works with former LIMA version

        @return: dict with interpreted metadata
        """
        metadata = {}

        if "entry" in self.input:
            self.entry = self.input_nxs.get_entry(self.input["entry"])
        else:
            self.entry = self.input_nxs.get_entries()[0]  # take the last entry
        instrument = self.input_nxs.get_class(self.entry, class_type="NXinstrument")
        if len(instrument) == 1:
            instrument = instrument[0]
        else:
            self.log_error("Expected ONE instrument is expected in entry, got %s in %s %s" %
                           (len(instrument), self.image_file, self.entry))
        detector_grp = self.input_nxs.get_class(instrument, class_type="NXdetector")
        if len(detector_grp) == 1:
            detector_grp = detector_grp[0]
        elif len(detector_grp) == 0 and "detector" in instrument:
            detector_grp = instrument["detector"]
        else:
            self.log_error("Expected ONE detector is expected in experiment, got %s in %s %s %s" %
                           (len(detector_grp), self.input_nxs, self.image_file, instrument))
        self.images_ds = detector_grp.get("data")
        if isinstance(self.images_ds, h5py.Group):
            if self.images_ds.attrs.get("NX_class") in ("NXdata", b"NXdata"):
                # this is an NXdata not a dataset: use the @signal
                self.images_ds = self.images_ds.get(self.images_ds.attrs.get("signal"))
        self.in_shape = self.images_ds.shape
        if "detector_information" in detector_grp:
            detector_information = detector_grp["detector_information"]
            if "name" in detector_information:
                metadata["DetectorName"] = ensure_str(detector_information["name"])
        # now read an interpret all static metadata.
        # metadata are on the collection side not instrument

        collections = self.input_nxs.get_class(self.entry, class_type="NXcollection")
        if len(collections) == 1:
            collection = collections[0]
        else:
            if len(collections) >= 1:
                collection = collections[0]
            else:
                self.log_error("Expected ONE collections is expected in entry, got %s in %s %s" %
                               (len(collections), self.image_file, self.entry))

        detector_grps = self.input_nxs.get_class(collection, class_type="NXdetector")
        if (len(detector_grps) == 0) and ("detector" in collection):
            detector_grp = collection["detector"]
        elif len(detector_grps) > 0:
            detector_grp = detector_grps[0]
        else:
            return {}
        if "parameters" in detector_grp:
            parameters_grp = detector_grp["parameters"]
        elif "data" in detector_grp:
            nxdata = detector_grp["data"]
            if "header" in nxdata and isinstance(nxdata["header"], h5py.Group):
                parameters_grp = nxdata["header"]
        else:
            return {}
        for key, value in parameters_grp.items():
            base = key.split("_")[0]
            conv = self.KEY_CONV.get(base, str)
            metadata[key] = conv(value[()])
        return metadata

    def parse_image_file_lima(self):
        """LIMA version working with LIMA-1.9.3

        @return: dict with interpreted metadata
        """
        metadata = {}

        if "entry" in self.input:
            self.entry = self.input_nxs.get_entry(self.input["entry"])
        else:
            self.entry = self.input_nxs.get_entries()[0]  # take the last entry
        instrument = self.input_nxs.get_class(self.entry, class_type="NXinstrument")
        if len(instrument) == 1:
            instrument = instrument[0]
        else:
            self.log_error("Expected ONE instrument is expected in entry, got %s in %s %s" %
                           (len(instrument), self.image_file, self.entry))
        detector_grp = self.input_nxs.get_class(instrument, class_type="NXdetector")
        if len(detector_grp) == 1:
            detector_grp = detector_grp[0]
        elif len(detector_grp) == 0 and "detector" in instrument:
            detector_grp = instrument["detector"]
        else:
            self.log_error("Expected ONE detector is expected in experiment, got %s in %s %s %s" %
                           (len(detector_grp), self.input_nxs, self.image_file, instrument))
        self.images_ds = detector_grp.get("data")
        if isinstance(self.images_ds, h5py.Group):
            if self.images_ds.attrs.get("NX_class") == "NXdata":
                # this is an NXdata not a dataset: use the @signal
                self.images_ds = self.images_ds.get(self.images_ds.attrs.get("signal"))
        self.in_shape = self.images_ds.shape
        if "detector_information" in detector_grp:
            detector_information = detector_grp["detector_information"]
            if "model" in detector_information:
                metadata["DetectorName"] = ensure_str(detector_information["model"][()])
        # now read an interpret all static metadata.
        # metadata are on the entry/instrument/detector/collection
        collections = self.input_nxs.get_class(detector_grp, class_type="NXcollection")
        headers = [grp for grp in collections if posixpath.basename(grp.name) == "header"]
        if len(headers) != 1:
            self.log_error("Expected a 'header' NXcollection in detector")
        parameters_grp = headers[0]
        for key, value in parameters_grp.items():
            base = key.split("_")[0]
            conv = self.KEY_CONV.get(base, ensure_str)
            metadata[key] = conv(value[()])
        return metadata

    def create_hdf5(self):
        """
        Create one HDF5 file per output
        Also initialize all workers
        """
        basename = os.path.splitext(os.path.basename(self.image_file))[0]
        if basename.endswith("_raw"):
            basename = basename[:-4]
        isotime = get_isotime()
        detector_grp = self.input_nxs.find_detector(all=True)
        detector_name = "undefined"
        for grp in detector_grp:
            if "detector_information/name" in grp:
                detector_name = grp["detector_information/name"][()]
        md_entry = self.metadata_nxs.get_entries()[0]
        instruments = self.metadata_nxs.get_class(md_entry, "NXinstrument")
        if instruments:
            collections = self.metadata_nxs.get_class(instruments[0], "NXcollection")
            to_copy = collections + detector_grp
        else:
            to_copy = detector_grp

        for ext in self.to_save:
            if ext == "raw":
                continue
            outfile = os.path.join(self.dest, "%s_%s.h5" % (basename, ext))
            self.output_hdf5[ext] = outfile
            try:
                nxs = Nexus(outfile, mode="a", creator="dahu")
            except IOError as error:
                self.log_warning("invalid HDF5 file %s: remove and re-create!\n%s" % (outfile, error))
                os.unlink(outfile)
                nxs = Nexus(outfile, mode="w", creator="dahu")
            entry = nxs.new_entry("entry",
                                  program_name=self.input.get("plugin_name", fully_qualified_name(self.__class__)),
                                  title=self.image_file + ":" + self.images_ds.name)
            entry["program_name"].attrs["version"] = __version__

            # Configuration
            config_grp = nxs.new_class(entry, "configuration", "NXnote")
            config_grp["type"] = "text/json"
            config_grp["data"] = json.dumps(self.input, indent=2, separators=(",\r\n", ": "))

            entry["detector_name"] = ensure_str(detector_name)

            nxprocess = nxs.new_class(entry, "PyFAI", class_type="NXprocess")
            nxprocess["program"] = "PyFAI"
            nxprocess["version"] = ensure_str(pyFAI.version)
            nxprocess["date"] = isotime
            nxprocess["processing_type"] = ensure_str(ext)
            nxdata = nxs.new_class(nxprocess, "result_" + ext, class_type="NXdata")
            entry.attrs["default"] = nxdata.name
            metadata_grp = nxprocess.require_group("parameters")

            for key, val in self.metadata.items():
                if isinstance(val, StringTypes):
                    metadata_grp[key] = ensure_str(val)
                else:
                    metadata_grp[key] = val

            # copy metadata from other files:
            for grp in to_copy:
                grp_name = posixpath.split(grp.name)[-1]
                if grp_name not in nxdata:
                    toplevel = nxprocess.require_group(grp_name)
                    for k, v in grp.attrs.items():
                        toplevel.attrs[k] = v
                else:
                    toplevel = nxprocess[grp_name]

                def grpdeepcopy(name, obj):
                    nxs.deep_copy(name, obj, toplevel=toplevel, excluded=["data"])

                grp.visititems(grpdeepcopy)

            shape = self.in_shape[:]
            if self.npt1_rad is None and "npt1_rad" in self.input:
                self.npt1_rad = int(self.input["npt1_rad"])
            else:
                qmax = self.ai.qArray(self.in_shape[-2:]).max()
                dqmin = self.ai.deltaQ(self.in_shape[-2:]).min() * 2.0
                self.npt1_rad = int(qmax / dqmin)

            if ext == "azim":
                if "npt2_rad" in self.input:
                    self.npt2_rad = int(self.input["npt2_rad"])
                else:
                    qmax = self.ai.qArray(self.in_shape[-2:]).max()
                    dqmin = self.ai.deltaQ(self.in_shape[-2:]).min() * 2.0
                    self.npt2_rad = int(qmax / dqmin)

                if "npt2_azim" in self.input:
                    self.npt2_azim = int(self.input["npt2_azim"])
                else:
                    chi = self.ai.chiArray(self.in_shape[-2:])
                    self.npt2_azim = int(numpy.degrees(chi.max() - chi.min()))
                shape = (self.in_shape[0], self.npt2_azim, self.npt2_rad)

                ai = self.ai.__copy__()
                ai.engines.update(ai.engines)  # copy engines as well

                worker = Worker(ai, self.in_shape[-2:], (self.npt2_azim, self.npt2_rad), self.unit)
                if self.flat is not None:
                    worker.ai.set_flatfield(self.flat)
                if self.dark is not None:
                    worker.ai.set_darkcurrent(self.dark)
                worker.output = "numpy"
                # if self.in_shape[0] < 5:
                #     worker.method = ("full", "histogram", "cython")  # "splitbbox"
                # else:
                worker.method = ("full", "csr", "opencl")  # "ocl_csr_gpu"

                if self.correct_solid_angle:
                    worker.set_normalization_factor(self.ai.pixel1 * self.ai.pixel2 / self.ai.dist / self.ai.dist)
                else:
                    worker.set_normalization_factor(1.0)
                    worker.correct_solid_angle = self.correct_solid_angle
                # self.log_warning("Normalization factor: %s" % worker.normalization_factor)

                worker.dummy = self.dummy
                worker.delta_dummy = self.delta_dummy
                if self.input.get("do_polarization"):
                    worker.polarization_factor = self.input.get("polarization_factor")
                worker.update_processor()
                self.workers[ext] = worker
            elif ext.startswith("ave"):
                if "_" in ext:
                    unit = ext.split("_", 1)[1]
                    npt1_rad = self.input.get("npt1_rad_" + unit, self.npt1_rad)
                    ai = self.ai.__copy__()
                else:
                    unit = self.unit
                    npt1_rad = self.npt1_rad
                    ai = self.ai
                shape = (self.in_shape[0], npt1_rad)
                worker = Worker(ai, self.in_shape[-2:], (1, npt1_rad), unit=unit)
                worker.output = "numpy"
                # if self.in_shape[0] < 5:
                #     worker.method = ("full", "histogram", "cython")  # "splitbbox"
                # else:
                worker.method = ("full", "csr", "opencl")  # "ocl_csr_gpu"
                if self.correct_solid_angle:
                    worker.set_normalization_factor(self.ai.pixel1 * self.ai.pixel2 / self.ai.dist / self.ai.dist)
                else:
                    worker.set_normalization_factor(1.0)
                    worker.correct_solid_angle = self.correct_solid_angle
                worker.dummy = self.dummy
                worker.delta_dummy = self.delta_dummy
                if self.input.get("do_polarization"):
                    worker.polarization_factor = True
                worker.update_processor()
                self.workers[ext] = worker
            elif ext == "sub":
                worker = PixelwiseWorker(dark=self.dark,
                                         dummy=self.dummy,
                                         delta_dummy=self.delta_dummy,
                                         )
                self.workers[ext] = worker
            elif ext == "flat":
                worker = PixelwiseWorker(dark=self.dark,
                                         flat=self.flat,
                                         dummy=self.dummy, delta_dummy=self.delta_dummy,
                                         )
                self.workers[ext] = worker
            elif ext == "solid":
                worker = PixelwiseWorker(dark=self.dark,
                                         flat=self.flat,
                                         solidangle=self.get_solid_angle(),
                                         dummy=self.dummy,
                                         delta_dummy=self.delta_dummy,
                                         polarization=self.polarization,
                                         )
                self.workers[ext] = worker
            elif ext == "dist":
                worker = DistortionWorker(dark=self.dark,
                                          flat=self.flat,
                                          solidangle=self.get_solid_angle(),
                                          dummy=self.dummy,
                                          delta_dummy=self.delta_dummy,
                                          polarization=self.polarization,
                                          detector=self.ai.detector,
                                          )
                self.workers[ext] = worker
                if self.distortion is None:
                    self.distortion = worker.distortion
                    self.cache_dis = str(self.ai.detector)
                    if self.cache_dis in self.cache:
                        self.distortion.lut = self.cache[self.cache_dis]
                    else:
                        self.distortion.calc_LUT()
                        self.cache[self.cache_dis] = self.distortion.lut
                else:
                    worker.distortion = self.distortion

            elif ext == "norm":
                worker = DistortionWorker(dark=self.dark,
                                          flat=self.flat,
                                          solidangle=self.get_solid_angle(),
                                          dummy=self.dummy,
                                          delta_dummy=self.delta_dummy,
                                          polarization=self.polarization,
                                          detector=self.ai.detector,
                                          mask=self.distortion_mask,
                                          device="gpu",
                                          method="csr"
                                          )
                self.workers[ext] = worker
                if self.distortion is None and worker.distortion is not None:
                    self.distortion = worker.distortion
                    self.cache_dis = str(self.ai.detector)
                    if self.cache_dis in self.cache:
                        self.distortion.lut = self.cache[self.cache_dis]
                    else:
                        self.distortion.calc_LUT()
                        self.cache[self.cache_dis] = self.distortion.lut
                else:
                    worker.distortion = self.distortion
            else:
                self.log_warning("unknown treatment %s" % ext)

            if (len(shape) >= 3):
                compression = {k: v for k, v in COMPRESSION.items()}
            else:
                compression = {}
            output_ds = nxdata.create_dataset("data",
                                              shape,
                                              dtype=numpy.float32,
                                              chunks=(1,) + shape[1:],
                                              maxshape=(None,) + shape[1:],
                                              **compression)
            nxdata.attrs["signal"] = "data"
            # output_ds.attrs["signal"] = "1"
            entry.attrs["default"] = nxdata.name
            if self.variance_formula is not None:
                error_ds = nxdata.create_dataset("data_errors", shape,
                                                 dtype=numpy.float32,
                                                 chunks=(1,) + shape[1:],
                                                 maxshape=(None,) + shape[1:],
                                                 **compression)
                # nxdata.attrs["uncertainties"] = "errors"
                self.output_ds[ext + "_err"] = error_ds
            if self.t is not None:
                nxdata["t"] = self.t
                nxdata["t"].attrs["interpretation"] = "scalar"
                nxdata["t"].attrs["unit"] = "s"

            if ext == "azim":
                nxdata.attrs["axes"] = [".", "chi", "q"]
                output_ds.attrs["interpretation"] = "image"
                if self.variance_formula is not None:
                    error_ds.attrs["interpretation"] = "image"

            elif ext == "ave":
                nxdata.attrs["axes"] = [".", "q"]
                output_ds.attrs["interpretation"] = "spectrum"
                if self.variance_formula is not None:
                    error_ds.attrs["interpretation"] = "spectrum"
            else:
                output_ds.attrs["interpretation"] = "image"
                if self.variance_formula is not None:
                    error_ds.attrs["interpretation"] = "image"

            self.output_ds[ext] = output_ds

    def process_images(self):
        """
        Here we process images....
        """
        if self.correct_I1:
            I1s = self.I1
        else:
            I1s = numpy.ones_like(self.I1)
        for i, I1 in enumerate(I1s):
            data = self.images_ds[i]
            if self.variance_formula:
                variance = self.variance_function(data, 0 if self.dark is None else self.dark)
            else:
                variance = None
            I1_corrected = I1 / self.scaling_factor
            for meth in self.to_save:
                if meth in ["raw", "dark"]:
                    continue
                res = err = None

                ds = self.output_ds[meth]
                if variance is not None:
                    err_ds = self.output_ds[meth + "_err"]

                if meth in ("sub", "flat", "dist", "cor"):
                    res = self.workers[meth].process(data, variance=variance)
                    if variance is not None:
                        res, err = res
                elif meth == "norm":
                    res = self.workers[meth].process(data, variance=variance,
                                                     normalization_factor=I1_corrected)
                    if variance is not None:
                        res, err = res[0], res[1]

                elif meth == "azim":
                    res = self.workers[meth].process(data, variance=variance,
                                                     normalization_factor=I1_corrected)
                    if (variance is not None):
                        if len(res) == 2:
                            res, err = res
                        else:
                            err = numpy.zeros_like(res)

                    if i == 0:
                        if "q" not in ds.parent:
                            ds.parent["q"] = numpy.ascontiguousarray(self.workers[meth].radial, dtype=numpy.float32)
                            ds.parent["q"].attrs["unit"] = self.unit
                            ds.parent["q"].attrs["axis"] = "3"
                            ds.parent["q"].attrs["interpretation"] = "scalar"
                        if "chi" not in ds.parent:
                            ds.parent["chi"] = numpy.ascontiguousarray(self.workers[meth].azimuthal, dtype=numpy.float32)
                            ds.parent["chi"].attrs["unit"] = "deg"
                            ds.parent["chi"].attrs["axis"] = "2"
                            ds.parent["chi"].attrs["interpretation"] = "scalar"
                elif meth.startswith("ave"):
                    res = self.workers[meth].process(data, variance=variance,
                                                     normalization_factor=I1_corrected)
                    # TODO: add other units
                    if i == 0 and "q" not in ds.parent:
                        if "log(1+q.nm)" in meth:
                            q = numpy.exp(self.workers[meth].radial) - 1.0
                        elif "log(1+q.A)" in meth:
                            q = (numpy.exp(self.workers[meth].radial) - 1.0) * 0.1
                        elif "log(q.nm)" in meth:
                            q = numpy.exp(self.workers[meth].radial)
                        elif "log10(q.m)" in meth:
                            q = 10 ** (self.workers[meth].radial) * 1e-9
                        else:
                            q = self.workers[meth].radial
                        ds.parent["q"] = numpy.ascontiguousarray(q, dtype=numpy.float32)
                        ds.parent["q"].attrs["unit"] = self.unit
                        ds.parent["q"].attrs["axis"] = "2"
                        ds.parent["q"].attrs["interpretation"] = "scalar"

                    if (variance is not None) and (res.shape[1] == 3):
                        err = res[:, 2]
                    else:
                        err = None
                    res = res[:, 1]
                else:
                    self.log_warning("Unknown/supported method ... %s, skipping" % (meth))
                    continue
                ds[i] = res
                if err is not None:
                    err_ds[i] = err

    def read_data(self, filename):
        """read dark data from a file

        @param filename: HDF5 file containing dark frames
        @return: numpy array with dark
        """
        if not os.path.exists(filename):
            self.log_error("Unable to read dark data from filename %s" % filename,
                           do_raise=True)
        with Nexus(filename, "r") as nxs:
            for entry in nxs.get_entries():
                for instrument in nxs.get_class(entry, "NXinstrument"):
                    for detector in nxs.get_class(instrument, "NXdetector"):
                        data = detector.get("data")
                        if isinstance(data, h5py.Group):
                            if data.attrs.get("NX_class") in ("NXdata", b"NXdata"):
                                # this is an NXdata not a dataset: use the @signal
                                data = data.get(data.attrs.get("signal"))
                        return numpy.array(data)

    def get_solid_angle(self):
        """ calculate the solid angle if needed and return it
        """
        if (self.absolute_solid_angle is None) and self.correct_solid_angle:
            self.absolute_solid_angle = self.ai.solidAngleArray(self.in_shape[-2:], absolute=True)
        return self.absolute_solid_angle

    def teardown(self):
        """Method for cleaning up stuff
        """
        # Finally update the cache:
        to_cache_array = {}
        to_cache_engine = {}
        for worker in self.workers.values():
            if "ai" in dir(worker):
                to_cache_array.update(worker.ai._cached_array)
                to_cache_engine.update(worker.ai.engines)

        cache_value = self.cache.get(self.cache_ai, CacheValues({}))
        cache_value.update(array=to_cache_array)
        cache_value.update(engine=to_cache_engine)

        # close also the source
        to_close = {}
        self.output_ds["source"] = self.images_ds
        for key, ds in self.output_ds.items():
            if not bool(ds):
                # the dataset is None when the file has been closed
                continue
            try:
                hdf5_file = ds.file
                filename = hdf5_file.filename
            except (RuntimeError, ValueError) as err:
                self.log_warning(f"Unable to retrieve filename of dataset {key}: {err}")
            else:
                to_close[filename] = hdf5_file
        for filename, hdf5_file in to_close.items():
            try:
                hdf5_file.close()
            except (RuntimeError, ValueError) as err:
                self.log_warning(f"Issue in closing file {filename} {type(err)}: {err}")

        self.ai = None
        self.polarization = None
        self.output["files"] = self.output_hdf5
        Plugin.teardown(self)
