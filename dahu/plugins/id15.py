#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Data Analysis plugin tailored for ID15

 
"""


from __future__ import with_statement, print_function


__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "23/03/2018"
__status__ = "development"
version = "0.4.0"

import os
import numpy
import logging
import copy
from dahu import version as dahu_version
from dahu.plugin import Plugin
from dahu.factory import register
from dahu.cache import DataCache
from dahu.utils import get_isotime
from threading import Thread, Event
import json
logger = logging.getLogger("id15")

try:
    from queue import Queue
except:
    from Queue import Queue
try:
    import pyFAI
except ImportError:
    logger.error("Failed to import PyFAI: download and install it from pypi")
try:
    import fabio
except ImportError:
    logger.error("Failed to import Fabio: download and install it from pypi")


class Reader(Thread):
    """Reader with input and output queue 
    """
    def __init__(self, queue_in, queue_out, quit_event):
        """Constructor of the class
        
        :param queue_in: input queue with (index, filename to read) as input
        :param queue_out: output queue with (index, FabioImage) as output
        :param quit_event: the event which tells the thread to end 
        """
        Thread.__init__(self)
        self._queue_in = queue_in
        self._queue_out = queue_out
        self._quit_event = quit_event

    def run(self):
        while not self._quit_event.is_set():
            plop = self._queue_in.get()
            if plop is None:
                break
            idx, fname = plop
            fimg = fabio.cbfimage.CbfImage()
            try:
                fimg = fimg.read(fname, check_MD5=False)
            except Exception as err:
                logger.error(err)
            self._queue_out.put((idx, fimg.data))
            self._queue_in.task_done()
            fimg = None

    @classmethod
    def build_pool(cls, args, size=1):
        """Create a pool of worker of a given size. 
        
        :param args: arguments to be passed to each of the worker
        :param size: size of the pool
        :return: a list of worker 
        """
        workers = []
        for _ in range(size):
            w = cls(*args)
            w.start()
            workers.append(w)
        return workers


class RawSaver(Thread):
    def __init__(self, queue, quit_event, dataset):
        """Constructor of the class
        
        :param queue: input queue with (index, raw_data) as input
        :param queue_out: output queue with (index, FabioImage) as output
        :param quit_event: the event which tells the thread to end
        :param dataset: the hdf5 dataset where to write 
        """
        Thread.__init__(self)
        self._queue = queue
        self._quit_event = quit_event
        self._dataset = dataset

    def run(self):
        while not self._quit_event.is_set():
            plop = self._queue.get()
            if plop is None:
                break
            idx, data = plop
            try:
                self._dataset[idx] = data
            except Exception as err:
                logger.error(err)
            self._queue.task_done()


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
    _ais = DataCache(10)  # key: poni-filename, value: ai
    pool_size = 6  # Default number of reader: needs to be tuned.

    def __init__(self):
        """
        """
        Plugin.__init__(self)
        self.queue_in = Queue()
        self.queue_out = Queue()
        self.queue_saver = None
        self.quit_event = Event()
        self.pool = []
        self.raw_saver = None

        self.ai = None  # this is the azimuthal integrator to use
        self.npt = 2000
        self.npt_azim = 256
        self.input_files = []
        self.method = "full_ocl_csr"
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
        ai = pyFAI.load(poni_file)
        stored = self._ais.get(poni_file, ai)
        if stored is ai:
            self.ai = stored
        else:
            self.ai = copy.deepcopy(stored)

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
            dataset = self.prepare_raw_hdf5(self.raw_compression)
            self.queue_saver = Queue()
            self.raw_saver = RawSaver(self.queue_saver, self.quit_event, dataset)
            self.raw_saver.start()
        # create the pool of workers
        self.pool = Reader.build_pool((self.queue_in, self.queue_out, self.quit_event), self.pool_size)

    def process(self):
        Plugin.process(self)
        logger.debug("IntegrateManyFrames.process")
        for idx, fname in enumerate(self.input_files):
            self.queue_in.put((idx, fname))
        if self.integration_method == "integrate2d":
            res = numpy.zeros((len(self.input_files), self.npt_azim, self.npt), dtype=numpy.float32)  # numpy array for storing data
        else:
            res = numpy.zeros((len(self.input_files), self.npt), dtype=numpy.float32)  # numpy array for storing data
        sigma = None
        if self.error_model:
            if self.integration_method == "integrate2d":
                sigma = numpy.zeros((len(self.input_files), self.npt_azim, self.npt), dtype=numpy.float32)  # numpy array for storing data
            else:
                sigma = numpy.zeros((len(self.input_files), self.npt), dtype=numpy.float32)  # numpy array for storing data

        method = self.ai.__getattribute__(self.integration_method)
        common_param = {"method": self.method,
                        "unit": self.unit,
                        "safe": False,
                        "dummy": self.dummy,
                        "delta_dummy": self.delta_dummy,
                        "error_model": self.error_model,
                        "mask": self.mask,
                        "polarization_factor": self.polarization_factor,
                        "normalization_factor": self.norm,
                        "correctSolidAngle": self.do_SA}
        if self.integration_method in ("integrate1d", "integrate_radial"):
            common_param["ntp"] = self.npt
        else:
            common_param["ntp_rad"] = self.npt
            common_param["ntp_azim"] = self.npt_azim
        if self.integration_method == "sigma_clip":
            common_param["thres"] = self.sigma_clip_thresold,
            common_param["max_iter"] = self.sigma_clip_max_iter
        if self.integration_method == "medfilt1d":
            common_param["percentile"] = self.medfilt1d_percentile

        for i in self.input_files:
            logger.debug("process %s", i)
            idx, data = self.queue_out.get()
            if data is None:
                self.log_error("Failed reading file: %s" % self.input_files[idx],
                               do_raise=False)
                continue
            if self.save_raw:
                self.queue_saver.put((idx, data))

            out = method(data, **common_param)
            res[idx] = out.intensity
            if self.error_model:
                sigma[idx] = out.sigma
            self.queue_out.task_done()

        self.queue_in.join()
        self.queue_out.join()
        if self.queue_saver is not None:
            self.queue_saver.join()

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
        coll = nxs.new_class(subentry, "process_integrate1d", class_type="NXdata")
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

        # now clean up threads, empty pool of workers
        self.quit_event.set()
        for _ in self.pool:
            self.queue_in.put(None)
        if self.queue_saver is not None:
            self.queue_saver.put(None)
            self.queue_saver = None
            self.raw_saver = None
        self.pool = None
        self.queue_in = None
        self.queue_out = None
        self.quit_event = None
