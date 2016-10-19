#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import with_statement, print_function

"""Data Analysis plugin tailored for ID15

 
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
from ..utils import get_isotime
from threading import Semaphore, Thread, Event
import logging
import copy
logger = logging.getLogger("id15")
import json
try:
    from queue import Queue
except:
    from Queue import Queue
try:
    import pyFAI
    from pyFAI.worker import make_ai
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
        """onstructor of the class
        
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
    
    Input parameters:
    :param poni_file: configuration of the geometry
    :param input_files:
    :param output_file:
    monitor_values
    
    Typical JSON file:
    {"poni_file": "/tmp/example.poni",
     "input_files": ["/tmp/file1.cbf", "/tmp/file2.cbf"],
     "monitor_values": [1, 1.1],
     "npt": 2000,
     "unit": "2th_deg",
     "output_file": "/path/to/dest.h5",
    }
    """
    _ais = DataCache(10)  # key: str(a), value= ai
    pool_size = 6 #Default number of reader: needs to be tuned. 

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
        self.ntp = 2000
        self.input_files = []
        self.method = "full_ocl_csr"
        self.unit = "q_nm^-1"
        self.output_file = None
        self.mask = ""
        self.wavelength = None
        self.polarization_factor = None
        self.do_SA = False
        self.norm = 1e12
        
    def setup(self, kwargs):
        logger.debug("Integrate.setup")
        Plugin.setup(self, kwargs)

        if "input_files" not in self.input:
            self.log_error("input_files not in input, save in input directory")
        self.input_files = 

        if "output_file" not in self.input:
            self.log_error("output_file not in input, save in input directory")
            self.output_file = "output.h5"
            # TODO: set path properly
        else:
            self.output_file = os.path.abspath(self.input["output_file"])
        if not self.output_file.endswith(".h5"):
            self.output_file = self.output_file + ".h5"

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
        for idx, fname in enumerate(self.input_files):
            self.queue_in.put((idx, fname))
        res = numpy.zeros((len(self.input_files), self.npt))  # numpy array for storing data

        for _ in self.input_files:
            idx, fimg = self.queue_out.get()
            if fimg.data is None:
                self.log_error("Failed reading file: %s" % self.input_files[idx],
                               do_raise=False)
                continue
            q, i = self.ai.integrate1d(fimg.data, self.npt, method="full_csr_ocl", unit=self.unit, safe=False)
            res[idx, :] = i
            self.queue_out.task_done()

        self.queue_in.join()
        self.queue_out.join()

        # now clean up threads
        self.quit_event.set()
        for _ in self.pool:
            self.queue_in.put(None)
        return res

#
#             fimg = fabio.open(fname)
#             if self.wavelength is not None:
#                 monitor = self.getMon(fimg.header, self.wavelength) / self.norm
#             else:
#                 monitor = 1.0
#             self.ai.integrate1d(fimg.data, npt=self.npt, method=self.method,
#                                 safe=False,
#                                 filename=destination,
#                                 normalization_factor=monitor,
#                                 unit=self.unit,
#                                 dummy=self.dummy,
#                                 delta_dummy=self.delta_dummy,
#                                 polarization_factor=self.polarization_factor,
#                                 correctSolidAngle=self.do_SA
#                                 )
#             self.output_files.append(destination)

    def save_result(self, q, I):
        isotime = numpy.string_(get_isotime())
        outfile = os.path.join(self.dest, "%s_%s.h5" % (basename, ext))
        self.output_hdf5[ext] = outfile
        try:
            nxs = pyFAI.io.Nexus(outfile, "a")
        except IOError as error:
            self.log_warning("invalid HDF5 file %s: remove and re-create!\n%s" % (outfile, error))
            os.unlink(outfile)
            nxs = pyFAI.io.Nexus(outfile)
        entry = nxs.new_entry("entry", program_name="dahu", title=self.image_file + ":" + self.images_ds.name)

        entry["program_name"].attrs["version"] = dahu.version
        entry["plugin_name"] = numpy.string_(".".join((os.path.splitext(os.path.basename(__file__))[0], self.__class__.__name__)))
        entry["plugin_name"].attrs["version"] = version
        entry["input"] = numpy.string_(json_config)
        entry["detector_name"] = numpy.string_(detector_name)

        subentry = nxs.new_class(entry, "PyFAI", class_type="NXprocess")
        subentry["program"] = numpy.string_("PyFAI")
        subentry["version"] = numpy.string_(pyFAI.version)
        subentry["date"] = isotime
        subentry["processing_type"] = numpy.string_(ext)
        coll = nxs.new_class(subentry, "process_" + ext, class_type="NXdata")
        metadata_grp = coll.require_group("parameters")

        for key, val in self.metadata.iteritems():
            if type(val) in [str, unicode]:
                metadata_grp[key] = numpy.string_(val)
            else:
                metadata_grp[key] = val

        # copy metadata from other files:
        for grp in to_copy:
            grp_name = posixpath.split(grp.name)[-1]
            if grp_name not in coll:
                toplevel = coll.require_group(grp_name)
                for k, v in grp.attrs.items():
                    toplevel.attrs[k] = v
            else:
                toplevel = coll[grp_name]

            def grpdeepcopy(name, obj):
                nxs.deep_copy(name, obj, toplevel=toplevel, excluded=["data"])

            grp.visititems(grpdeepcopy)

        shape = self.in_shape[:]

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

            ai = self.ai.__deepcopy__()
            worker = pyFAI.worker.Worker(ai, self.in_shape[-2:], (self.npt2_azim, self.npt2_rad), "q_nm^-1")
            if self.flat is not None:
                worker.ai.set_flatfield(self.flat)
            if self.dark is not None:
                worker.ai.set_darkcurrent(self.dark)
            worker.output = "numpy"
            if self.in_shape[0] < 5:
                worker.method = "splitbbox"
            else:
                worker.method = "ocl_csr_gpu"
            if self.correct_solid_angle:
                worker.set_normalization_factor(self.ai.pixel1 * self.ai.pixel2 / self.ai.dist / self.ai.dist / self.scaling_factor)
            else:
                worker.set_normalization_factor(1.0 / self.scaling_factor)
                worker.correct_solid_angle = self.correct_solid_angle
            self.log_warning("Normalization factor: %s" % worker.normalization_factor)

            worker.dummy = self.dummy
            worker.delta_dummy = self.delta_dummy
            self.workers[ext] = worker
        elif ext == "ave":
            if "npt1_rad" in self.input:
                self.npt1_rad = int(self.input["npt1_rad"])
            else:
                qmax = self.ai.qArray(self.in_shape[-2:]).max()
                dqmin = self.ai.deltaQ(self.in_shape[-2:]).min() * 2.0
                self.npt1_rad = int(qmax / dqmin)
            shape = (self.in_shape[0], self.npt1_rad)
            worker = pyFAI.worker.Worker(self.ai, self.in_shape[-2:], (1, self.npt1_rad), "q_nm^-1")
            worker.output = "numpy"
            if self.in_shape[0] < 5:
                worker.method = "splitbbox"
            else:
                worker.method = "ocl_csr_gpu"
            if self.correct_solid_angle:
                worker.set_normalization_factor(self.ai.pixel1 * self.ai.pixel2 / self.ai.dist / self.ai.dist / self.scaling_factor)
            else:
                worker.set_normalization_factor(1.0 / self.scaling_factor)
                worker.correct_solid_angle = self.correct_solid_angle
            worker.dummy = self.dummy
            worker.delta_dummy = self.delta_dummy
            self.workers[ext] = worker
        elif ext == "sub":
            worker = pyFAI.worker.PixelwiseWorker(dark=self.dark)
            self.workers[ext] = worker
        elif ext == "flat":
            worker = pyFAI.worker.PixelwiseWorker(dark=self.dark, flat=self.flat)
            self.workers[ext] = worker
        elif ext == "solid":
            worker = pyFAI.worker.PixelwiseWorker(dark=self.dark, flat=self.flat, solidangle=self.get_solid_angle())
            self.workers[ext] = worker
        elif ext == "dist":
            worker = pyFAI.worker.DistortionWorker(dark=self.dark, flat=self.flat, solidangle=self.get_solid_angle(),
                                                   detector=self.ai.detector)
            self.workers[ext] = worker
        elif ext == "norm":
            worker = pyFAI.worker.DistortionWorker(dark=self.dark, flat=self.flat, solidangle=self.get_solid_angle(),
                                                   detector=self.ai.detector)
            self.workers[ext] = worker
        else:
            self.log_warning("unknown treatment %s" % ext)
        output_ds = coll.create_dataset("data", shape, "float32",
                                        chunks=(1,) + shape[1:],
                                        maxshape=(None,) + shape[1:])
        if self.t is not None:
            coll["t"] = self.t
            coll["t"].attrs["axis"] = "1"
            coll["t"].attrs["interpretation"] = "scalar"
            coll["t"].attrs["unit"] = "s"

#             output_ds.attrs["NX_class"] = "NXdata" -> see group
        output_ds.attrs["signal"] = "1"
        if ext == "azim":
            output_ds.attrs["axes"] = ["t", "chi", "q"]
            output_ds.attrs["interpretation"] = "image"
        elif ext == "ave":
            output_ds.attrs["axes"] = ["t", "q"]
            output_ds.attrs["interpretation"] = "spectrum"
        elif ext in ("sub", "flat", "solid", "dist"):
            output_ds.attrs["axes"] = "t"
            output_ds.attrs["interpretation"] = "image"
        else:
            output_ds.attrs["interpretation"] = "image"
        self.output_ds[ext] = output_ds


    def teardown(self):
        Plugin.teardown(self)
        logger.debug("Integrate.teardown")
        # Create some output data
        self.output["output_files"] = self.output_files

