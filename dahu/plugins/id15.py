#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import with_statement, print_function

"""Data Analysis plugin tailored for ID15

 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "25/10/2016"
__status__ = "development"
version = "0.1.0"

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
            self._queue_out.put((idx, fimg))
            self._queue_in.task_done()


@register
class IntegrateManyFrames(Plugin):
    """This is the basic plugin of PyFAI for azimuthal integration
       
    Typical JSON file:
    {"poni_file": "/tmp/example.poni",
     "input_files": ["/tmp/file1.cbf", "/tmp/file2.cbf"],
     "monitor_values": [1, 1.1],
     "npt": 2000,
     "unit": "2th_deg",
     "output_file": "/path/to/dest.h5",
     "mask":
     "wavelength"
     "unit"
     "dummy"
     "delta_dummy"
     "do_polarziation"
     "polarization_factor"
     "do_SA"
     "norm"
    }
    """
    _ais = DataCache(10)  # key: poni-filename, value: ai
    pool_size = 6  # Default number of reader: needs to be tuned.

    @staticmethod
    def build_pool(worker, args, size=1):
        """Create a pool of worker of a given size. 
        
        :param worker: class of the worker (deriving  from Thread) 
        :param args: arguments to be passed to each of the worker
        :param size: size of the pool
        :return: a list of worker 
        """
        workers = []
        for _ in range(size):
            w = worker(*args)
            w.start()
            workers.append(w)
        return workers

    def __init__(self):
        """
        """
        Plugin.__init__(self)
        self.queue_in = Queue()
        self.queue_out = Queue()
        self.quit_event = Event()
        self.pool = []

        self.ai = None  # this is the azimuthal integrator to use
        self.npt = 2000
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

    def setup(self, kwargs=None):
        """Perform the setup of the job.
        mainly parsing of the kwargs. 
         
        :param kwargs: dict with parmaters.
        :return: None
        """
        logger.debug("IntegrateManyFrames.setup")
        Plugin.setup(self, kwargs)

        # create the pool of workers
        self.pool = self.build_pool(Reader, (self.queue_in, self.queue_out, self.quit_event), self.pool_size)

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

    def process(self):
        Plugin.process(self)
        logger.debug("IntegrateManyFrames.process")
        for idx, fname in enumerate(self.input_files):
            self.queue_in.put((idx, fname))
        res = numpy.zeros((len(self.input_files), self.npt), dtype=numpy.float32)  # numpy array for storing data
        sigma = None
        if self.error_model:
            sigma = numpy.zeros((len(self.input_files), self.npt), dtype=numpy.float32)  # numpy array for storing data
        for i in self.input_files:
            logger.debug("process %s", i)
            idx, fimg = self.queue_out.get()
            if fimg.data is None:
                self.log_error("Failed reading file: %s" % self.input_files[idx],
                               do_raise=False)
                continue
            out = self.ai.integrate1d(fimg.data, self.npt, method=self.method,
                                      unit=self.unit, safe=False,
                                      dummy=self.dummy, delta_dummy=self.delta_dummy,
                                      error_model=self.error_model,
                                      mask=self.mask,
                                      polarization_factor=self.polarization_factor,
                                      normalization_factor=self.norm,
                                      correctSolidAngle=self.do_SA)
            res[idx, :] = out.intensity
            if self.error_model:
                sigma[idx, :] = out.sigma
            self.queue_out.task_done()

        self.queue_in.join()
        self.queue_out.join()

        self.save_result(out, res, sigma)

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

        # Explicit handling of list of unicode names in needed !
#         input_grp = entry.require_group("input")
#         for key, value in self.input.items():
#             try:
#                 if isinstance(value, (str, unicode)):
#                     input_grp[key] = numpy.string_(value)
#                 elif __len__ in dir(value) and
#                 else:
#                     input_grp[key] = value
#             except TypeError:
#                 logger.warning("Conversion error for %s: %s", key, value)
        # If more speed is needed:

        entry["input"] = numpy.string_(json.dumps(self.input))

        subentry = nxs.new_class(entry, "PyFAI", class_type="NXprocess")
        subentry["program"] = numpy.string_("PyFAI")
        subentry["version"] = numpy.string_(pyFAI.version)
        subentry["date"] = isotime
        subentry["processing_type"] = numpy.string_("integrate1d")
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

        # now clean up threads, empty pool of workers
        self.quit_event.set()
        for _ in self.pool:
            self.queue_in.put(None)
        self.pool = None
        self.queue_in = None
        self.queue_out = None
        self.quit_event = None
