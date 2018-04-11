#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Data Analysis plugin tailored for ID15

 
"""

from __future__ import with_statement, print_function

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "11/04/2018"
__status__ = "development"
version = "0.4.0"

import os
import numpy
import logging
import copy
import json
import glob
from dahu import version as dahu_version
from dahu.plugin import Plugin
from dahu.factory import register
from dahu.utils import get_isotime

logger = logging.getLogger("id15v2")

try:
    import fabio
    import pyFAI
    from silx.opencl.codec.byte_offset import ByteOffset
except ImportError:
    logger.error("Failed to import PyFAI, fabio or silx: download and install it from pypi")


@register
class IntegrateManyFrames(Plugin):
    """This is the basic plugin of PyFAI for azimuthal integration
       
    Typical JSON file:
    {"poni_file": "/tmp/example.poni",
     "input_files": ["/tmp/file1.cbf", "/tmp/file2.cbf"],
     "monitor_values": [1, 1.1],
     "npt": 2000,
     "npt_azim": None, 
     "unit": "2th_deg",
     "output_file": "/path/to/dest.h5",
     "save_raw": "/path/to/hdf5/with/raw.h5",
     "delete_incoming": False,
     "do_SA":True,
     }
     Plus:
     "mask":
     "wavelength"
     "unit"
     "dummy"
     "delta_dummy"
     "do_polarziation"
     "polarization_factor"
     "do_SA"
     "norm"
     "raw_compression":  "bitshuffle",
     "integration_method": "integrate1d"
     "sigma_clip_thresold":3
     self.sigma_clip_max_iter: 5
    
    """

    def __init__(self):
        """
        """
        Plugin.__init__(self)
        self.ai = None  # this is the azimuthal integrator to use
        self.npt = 2000
        self.npt_azim = 256
        self.input_files = []
        self.method = "ocl_nosplit_csr_gpu"
        self.unit = "q_nm^-1"
        self.output_file = None
        self.mask = None
        self.wavelength = None
        self.polarization_factor = None
        self.do_SA = False
        self.dummy = None
        self.delta_dummy = None
        self.norm = 1
        self.error_model = None  # "poisson"
        self.save_raw = None
        self.raw_nxs = None
        self.raw_ds = None
        self.raw_compression = None
        self.integration_method = "integrate1d"
        self.sigma_clip_thresold = 3
        self.sigma_clip_max_iter = 5
        self.medfilt1d_percentile = (10, 90)

    def setup(self, kwargs=None):
        """Perform the setup of the job.
        mainly parsing of the kwargs. 
         
        :param kwargs: dict with parmaters.
        :return: None
        """
        logger.debug("IntegrateManyFrames.setup")
        Plugin.setup(self, kwargs)

        self.input_files = self.input.get("input_files")
        if not self.input_files:
            self.log_error("InputError: input_files not in input.")
        if not isinstance(self.input_files, list):
            self.input_files = glob.glob(self.input_files)
            self.input_files.sort()
        if "output_file" not in self.input:
            self.log_error("InputWarning: output_file not in input, save in input directory",
                           do_raise=False)
            self.output_file = os.path.join(os.path.dirname(self.input_files[0]), "output.h5")
        else:
            self.output_file = os.path.abspath(self.input["output_file"])
        if not self.output_file.endswith(".h5"):
            self.output_file = self.output_file + ".h5"

        poni_file = self.input.get("poni_file")
        if not poni_file:
            self.log_error("InputError: poni_file not in input.")
        self.ai = pyFAI.load(poni_file)
#         stored = self._ais.get(poni_file, ai)
#         if stored is ai:
#             self.ai = stored
#         else:
#             self.ai = copy.deepcopy(stored)

        self.npt = int(self.input.get("npt", self.npt))
        self.npt_azim = self.input.get("npt_azim", self.npt_azim)
        self.unit = self.input.get("unit", self.unit)
        self.wavelength = self.input.get("wavelength", self.wavelength)
        if os.path.exists(self.input.get("mask", "")):
            self.mask = fabio.open(self.input["mask"]).data
        self.dummy = self.input.get("dummy", self.dummy)
        self.delta_dummy = self.input.get("delta_dummy", self.delta_dummy)
        if self.input.get("do_polarziation"):
            self.polarization_factor = self.input.get("polarization_factor", self.polarization_factor)
        self.do_SA = self.input.get("do_SA", self.do_SA)
        self.norm = self.input.get("norm", self.norm)
        self.method = self.input.get("method", self.method)
        self.save_raw = self.input.get("save_raw", self.save_raw)
        self.integration_method = self.input.get("integration_method", self.integration_method)
        self.sigma_clip_thresold = self.input.get("sigma_clip_thresold", self.sigma_clip_thresold)
        self.sigma_clip_max_iter = self.input.get("sigma_clip_max_iter", self.sigma_clip_max_iter)
        self.medfilt1d_percentile = self.input.get("medfilt1d_percentile", self.medfilt1d_percentile)

        self.raw_compression = self.input.get("raw_compression", self.raw_compression)
        if self.save_raw:
            self.prepare_raw_hdf5(self.raw_compression)

    def process(self):
        Plugin.process(self)
        logger.debug("IntegrateManyFrames.process")
        if self.integration_method == "integrate2d":
            res = numpy.zeros((len(self.input_files), self.npt_azim, self.npt), dtype=numpy.float32)  # numpy array for storing data
        else:
            res = numpy.zeros((len(self.input_files), self.npt), dtype=numpy.float32)  # numpy array for storing data
        sigma = None
        if self.error_model or self.integration_method == "sigma_clip":
            if self.integration_method == "integrate2d":
                sigma = numpy.zeros((len(self.input_files), self.npt_azim, self.npt), dtype=numpy.float32)  # numpy array for storing data
            else:
                sigma = numpy.zeros((len(self.input_files), self.npt), dtype=numpy.float32)  # numpy array for storing data

        method = self.ai.__getattribute__(self.integration_method)
        common_param = {"method": self.method,
                        "unit": self.unit,
                        "dummy": self.dummy,
                        "delta_dummy": self.delta_dummy,
                        "mask": self.mask,
                        "polarization_factor": self.polarization_factor,
                        "normalization_factor": self.norm,
                        "correctSolidAngle": self.do_SA}
        if self.integration_method in ("integrate1d", "integrate_radial"):
            common_param["npt"] = self.npt
            common_param["error_model"] = self.error_model
            common_param["safe"] = False
        else:
            common_param["npt_rad"] = self.npt
            common_param["npt_azim"] = self.npt_azim
        if self.integration_method == "sigma_clip":
            common_param["thres"] = self.sigma_clip_thresold,
            common_param["max_iter"] = self.sigma_clip_max_iter
        if self.integration_method == "medfilt1d":
            common_param["percentile"] = self.medfilt1d_percentile

        # prepare some tools
        cbf = fabio.open(self.input_files[0])
        bo = ByteOffset(os.path.getsize(self.input_files[0]), cbf.data.size,
                        devicetype="gpu")
        shape = cbf.data.shape
        for idx, fname in enumerate(self.input_files):
            logger.debug("process %s: %s", idx, fname)
            if fname.endswith("cbf"):
                raw = cbf.read(fname, only_raw=True)
                data = bo(raw, as_float=False).get().reshape(shape)
            else:
                data = fabio.open(fname).data
            if data is None:
                self.log_error("Failed reading file: %s" % self.input_files[idx],
                               do_raise=False)
                continue
            if self.save_raw:
                self.raw_ds[idx] = data

            out = method(data, **common_param)
            res[idx] = out.intensity
            if self.error_model or self.integration_method == "sigma_clip":
                sigma[idx] = out.sigma

        self.save_result(out, res, sigma)
        if self.input.get("delete_incoming"):
            for fname in self.input_files:
                try:
                    os.unlink(fname)
                except IOError as err:
                    self.log_warning(err)

    def prepare_raw_hdf5(self, filter_=None):
        """Prepare an HDF5 output file for saving raw data
        :param filter_: name of the compression filter 
        """
        kwfilter = {}
        if filter_ == "gzip":
            kwfilter = {"compression": "gzip", "shuffle": True}
        elif filter_ == "lz4":
            kwfilter = {"compression": 32004, "shuffle": True}
        elif filter_ == "bitshuffle":
            kwfilter = {"compression": 32008, "compression_opts": (0, 2)}  # enforce lz4 compression

        first_image = self.input_files[0]
        fimg = fabio.open(first_image)
        shape = fimg.data.shape
        stack_shape = (len(self.input_files),) + shape
        first_frame_timestamp = os.stat(first_image).st_ctime

        try:
            self.raw_nxs = pyFAI.io.Nexus(self.save_raw, "a")
        except IOError as error:
            self.log_warning("invalid HDF5 file %s: remove and re-create!\n%s" % (self.save_raw, error))
            os.unlink(self.save_raw)
            self.raw_nxs = pyFAI.io.Nexus(self.save_raw)
        entry = self.raw_nxs.new_entry("entry",
                                       program_name="dahu",
                                       title="ID15.raw_data",
                                       force_time=first_frame_timestamp,
                                       force_name=True)
        entry["program_name"].attrs["version"] = dahu_version
        entry["plugin_name"] = numpy.string_(".".join((os.path.splitext(os.path.basename(__file__))[0], self.__class__.__name__)))
        entry["plugin_name"].attrs["version"] = version
        coll = self.raw_nxs.new_class(entry, "data", class_type="NXdata")
        try:
            self.raw_ds = coll.require_dataset(name="data", shape=stack_shape,
                                               dtype=fimg.data.dtype,
                                               chunks=(1,) + shape,
                                               **kwfilter)
        except Exception as error:
            logger.error("Error in creating dataset, disabling compression:%s", error)
            self.raw_ds = coll.require_dataset(name="data", shape=stack_shape,
                                               dtype=fimg.data.dtype,
                                               chunks=(1,) + shape)
        return self.raw_ds

    def save_result(self, out, I, sigma=None):
        """Save the result of the work as a HDF5 file
        
        :param out: scattering result 
        :param I: Intensities as 2D array
        :param sigma: standard deviation of I as 2D array, is possible 
        """
        logger.debug("IntegrateManyFrames.save_result")
        isotime = numpy.string_(get_isotime())
        try:
            nxs = pyFAI.io.Nexus(self.output_file, "a")
        except IOError as error:
            self.log_warning("invalid HDF5 file %s: remove and re-create!\n%s" % (self.output_file, error))
            os.unlink(self.output_file)
            nxs = pyFAI.io.Nexus(self.output_file)
        entry = nxs.new_entry("entry", program_name="dahu", title="ID15.IntegrateManyFrames ")

        entry["program_name"].attrs["version"] = dahu_version
        entry["plugin_name"] = numpy.string_(".".join((os.path.splitext(os.path.basename(__file__))[0], self.__class__.__name__)))
        entry["plugin_name"].attrs["version"] = version

        entry["input"] = numpy.string_(json.dumps(self.input))
        entry["input"].attrs["format"] = 'json'
        subentry = nxs.new_class(entry, "PyFAI", class_type="NXprocess")
        subentry["program"] = numpy.string_("PyFAI")
        subentry["version"] = numpy.string_(pyFAI.version)
        subentry["date"] = isotime
        subentry["processing_type"] = numpy.string_(self.integration_method)
        coll = nxs.new_class(subentry, "process_%s" % self.integration_method,
                             class_type="NXdata")
        metadata_grp = coll.require_group("parameters")
        for key, value in self.ai.getPyFAI().items():
            metadata_grp[key] = numpy.string_(value)
        scale, unit = str(out.unit).split("_", 1)
        coll[scale] = out.radial.astype("float32")
        coll[scale].attrs["interpretation"] = "scalar"
        coll[scale].attrs["unit"] = unit

        coll["I"] = I.astype("float32")
        coll["I"].attrs["interpretation"] = "spectrum"
        coll["I"].attrs["axes"] = ["t", scale]
        coll["I"].attrs["signal"] = "1"

        if sigma is not None:
            coll["sigma"] = sigma.astype("float32")
            coll["sigma"].attrs["interpretation"] = "spectrum"
            coll["sigma"].attrs["axes"] = ["t", scale]
        nxs.close()

    def teardown(self):
        Plugin.teardown(self)
        logger.debug("IntegrateManyFrames.teardown")
        # Create some output data
        self.output["output_file"] = self.output_file
        if self.save_raw:
            self.raw_nxs.close()
            self.output["save_raw"] = self.save_raw
            self.raw_nxs = None
            self.raw_ds = None
