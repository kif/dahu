#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* Mesh mode: Rebuild the complete map and performs basic analysis on it.
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "28/04/2025"
__status__ = "development"
__version__ = "0.1.0"

import json
from dataclasses import dataclass, fields, asdict
import numpy
from dahu.plugin import Plugin
from .common import Sample, Ispyb, get_equivalent_frames, cmp_float, get_integrator, KeyCache, \
                    polarization_factor, method, Nexus, get_isotime, SAXS_STYLE, NORMAL_STYLE, \
                    create_nexus_sample

@dataclass(slots=True)
class Scan:
    """Class describing 2D mesh-scan"""
    fast_motor_name: str = "fast"
    fast_motor_start: float = 0.0
    fast_motor_stop: float = 1.0
    fast_motor_step: int = 0  # bliss-like steps, actually step+1 points
    slow_motor_name: str = "slow"
    slow_motor_start: float = 0.0
    slow_motor_stop: float = 1.0
    slow_motor_step: int = 0  # bliss-like steps, actually step+1 points
    backnforth: bool = False

    def __repr__(self):
        return json.dumps(self.as_dict(), indent=4)

    def as_dict(self):
        """Like asdict, without extra features:

        :return: dict which can be JSON-serialized 
        """
        dico = {}
        for key, value in asdict(self).items():
            dico[key] = value
        return dico

    @classmethod
    def from_dict(cls, dico):
        """Alternative constructor,
        :param dico: dict with the config
        :return: instance of the dataclass
        """
        to_init = dico
        self = cls(**to_init)
        return self

    def save(self, filename):
        """Dump the content of the dataclass as JSON file"""
        with open(filename, "w") as w:
            w.write(json.dumps(self.as_dict(), indent=2))




def input_from_master(master_file):
    """Convert a bliss masterfile containing a single NXentry with a 2D scan into a set of plugins to be launched

    :param master_file: path to a bliss masterfile
    :return: list of json dicts.
    """
    pass


class Mesh(Plugin):
    """ Rebuild the complete map and perform basic analysis on it.
    
        Typical JSON file:
    {
      "integrated_files": ["img_001.h5", "img_002.h5"],
      "output_file": "mesh.h5"
      "ispyb": {
        "url": "http://ispyb.esrf.fr:1234",
        "pyarch": "/data/pyarch/mx1234/sample", 
        "measurement_id": -1,
        "collection_id": -1
       },
       "nmf_components": 5, 
       "scan": {
            "fast_motor_name": "chipz",
            "fast_motor_start": -2.2,
            "fast_motor_stop": -3.2,
            "fast_motor_step": 3, 
            "slow_motor_name": "chipy",
            "slow_motor_start": -5.2,
            "slow_motor_stop": -12.2,
            "slow_motor_step": 7, 
            "backnforth": False,
            }
      "wait_for": [jobid_img001, jobid_img002],
      "plugin_name": "bm29.mesh"
    } 
    """
    
    def __init__(self):
        Plugin.__init__(self)
        self.input_files = []
        self.nxs = None
        self.output_file = None
        self.scan = Scan()
        self.juices = []
        self.to_pyarch = {}
        self.ispyb = None
        self._pid = 0

    def sequence_index(self):
        value = self._pid
        self._pid += 1
        return value

    def setup(self):
        Plugin.setup(self)

        for job_id in self.input.get("wait_for", []):
            self.wait_for(job_id)

        self.input_files = [os.path.abspath(i) for i in self.input.get("integrated_files", "")]
        self.output_file = self.input.get("output_file")
        if not self.output_file:
            dirname, basename = os.path.split(os.path.commonprefix(self.input_files) + "_mesh.h5")
            dirname = os.path.dirname(dirname)
            dirname = os.path.join(dirname, "mesh")
            self.output_file = os.path.join(dirname, basename)
            if not os.path.isdir(dirname):
                try:
                    os.makedirs(dirname)
                except Exception as err:
                    self.log_warning(f"Unable to create dir {dirname}. {type(err)}: {err}")

            self.log_warning(f"No output file provided, using {self.output_file}")

        #Manage gallery here
        dirname = os.path.dirname(self.output_file)
        gallery = os.path.join(dirname, "gallery")
        if not os.path.isdir(gallery):
            try:
                os.makedirs(gallery)
            except Exception as err:
                self.log_warning(f"Unable to create dir {gallery}. {type(err)}: {err}")
        ispydict = self.input.get("ispyb", {})
        ispydict["gallery"] = gallery
        self.ispyb = Ispyb._fromdict(ispydict)
        self.scan = Scan.from_dict(self.input.scan)

    def process(self):
        self.create_nexus()
        self.to_pyarch["hdf5_filename"] = self.output_file
        # self.to_pyarch["chunk_size"] = self.juices[0].Isum.size
        self.to_pyarch["id"] = os.path.commonprefix(self.input_files)
        self.to_pyarch["sample_name"] = self.juices[0].sample.name
        if not self.input.get("no_ispyb"):
            self.send_to_ispyb()
        # self.output["icat"] = 
        self.send_to_icat()

    def teardown(self):
        Plugin.teardown(self)
        logger.debug("HPLC.teardown")
        # export the output file location
        self.output["output_file"] = self.output_file
        if self.nxs is not None:
            self.nxs.close()
        self.to_pyarch = None
        self.ispyb = None

    def create_nexus(self):
        nxs = Nexus(self.output_file, mode="w")
        entry_grp = nxs.new_entry("entry", self.input.get("plugin_name", "dahu"),
                              title='BioSaxs Mesh experiment',
                              force_time=get_isotime())
        entry_grp["version"] = __version__
        nxs.h5.attrs["default"] = entry_grp.name

    # Configuration
        cfg_grp = nxs.new_class(entry_grp, "configuration", "NXnote")
        cfg_grp.create_dataset("data", data=json.dumps(self.input, indent=2, separators=(",\r\n", ":\t")))
        cfg_grp.create_dataset("format", data="text/json")

    # Process 0: Measurement group
        input_grp = nxs.new_class(entry_grp, "0_measurement", "NXcollection")
        input_grp["sequence_index"] = self.sequence_index()

        for idx, filename in enumerate(self.input_files):
            juice = self.read_nexus(filename)
            if juice is not None:
                rel_path = os.path.relpath(os.path.abspath(filename), os.path.dirname(os.path.abspath(self.output_file)))
                input_grp["LImA_%04i" % idx] = h5py.ExternalLink(rel_path, juice.h5path)
                self.juices.append(juice)

        q = self.juices[0].q
        unit = self.juices[0].unit
        radial_unit, unit_name = str(unit).split("_", 1)

        # Sample: outsourced !
        create_nexus_sample(nxs, entry_grp, self.juices[0].sample)

    # Process 1: Mesh
        mesh_grp = nxs.new_class(entry_grp, "1_mesh", "NXprocess")
        mesh_grp["sequence_index"] = self.sequence_index()
        nframes = max(i.idx.max() for i in self.juices) + 1
        nbin = q.size

        shape = (self.scan.slow_motor_step+1, self.scan.fast_motor_step+1, nbin)
        I = numpy.zeros(shape, dtype=numpy.float32)
        slow = numpy.linspace(self.scan.slow_motor_start,
                              self.scan.slow_motor_stop,
                              self.scan.slow_motor_step+1,
                              endpoint=True).astype("float32")
        fast = numpy.linspace(self.scan.fast_motor_start,
                              self.scan.fast_motor_stop,
                              self.scan.fast_motor_step+1,
                              endpoint=True).astype("float32")

        sigma = numpy.zeros(shape, dtype=numpy.float32)
        Isum = numpy.zeros(shape, dtype=numpy.float32)

        #TODO: order images accordingly:

        ids = numpy.arange(nframes)
        idx = numpy.concatenate([i.idx for i in self.juices])
        timestamps = self.to_pyarch["time"] = numpy.concatenate([i.timestamps for i in self.juices])
        I[idx] = numpy.vstack([i.I for i in self.juices])
        Isum[idx] = numpy.concatenate([i.Isum for i in self.juices])
        sigma[idx] = numpy.vstack([i.sigma for i in self.juices])

        hplc_data = nxs.new_class(mesh_grp, "hplc", "NXdata")
        hplc_data.attrs["title"] = "Chromatogram"
        sum_ds = hplc_data.create_dataset("sum", data=Isum, dtype=numpy.float32)
        sum_ds.attrs["interpretation"] = "spectrum"
        sum_ds.attrs["long_name"] = "Summed Intensity"
        frame_ds = hplc_data.create_dataset("frame_ids", data=ids, dtype=numpy.uint32)
        frame_ds.attrs["interpretation"] = "spectrum"
        frame_ds.attrs["long_name"] = "frame index"
        hplc_data.attrs["signal"] = "sum"
        hplc_data.attrs["axes"] = "frame_ids"
        mesh_grp.attrs["default"] = entry_grp.attrs["default"] = hplc_data.name
        time_ds = hplc_data.create_dataset("timestamps", data=timestamps, dtype=numpy.uint32)
        time_ds.attrs["interpretation"] = "spectrum"
        time_ds.attrs["long_name"] = "Time stamps (s)"

        integration_data = nxs.new_class(mesh_grp, "results", "NXdata")
        mesh_grp.attrs["title"] = str(self.juices[0].sample)

        int_ds = integration_data.create_dataset("I", data=numpy.ascontiguousarray(I, dtype=numpy.float32))
        std_ds = integration_data.create_dataset("errors", data=numpy.ascontiguousarray(sigma, dtype=numpy.float32))
        q_ds = integration_data.create_dataset("q", data=self.juices[0].q)
        q_ds.attrs["interpretation"] = "spectrum"
        q_ds.attrs["unit"] = unit_name
        q_ds.attrs["long_name"] = "Scattering vector q (nm⁻¹)"
        integration_data.attrs["signal"] = "I"
        integration_data.attrs["axes"] = [".", "q"]
        integration_data.attrs["SILX_style"] = SAXS_STYLE

        int_ds.attrs["interpretation"] = "spectrum"
        int_ds.attrs["units"] = "arbitrary"
        int_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"
        # int_ds.attrs["uncertainties"] = "errors" This does not work
        int_ds.attrs["scale"] = "log"
        std_ds.attrs["interpretation"] = "spectrum"

        save_zip(os.path.splitext(self.output_file)[0]+".zip", 
                 self.juices[0], I, sigma)
    
    @staticmethod
    def read_nexus(filename):
        "return some NexusJuice from a HDF5 file "
        with Nexus(filename, "r") as nxsr:
            entry_name = nxsr.h5.attrs["default"]
            entry_grp = nxsr.h5[entry_name]
            h5path = entry_grp.name
            nxdata_grp = nxsr.h5[entry_grp.attrs["default"]]
            # assert nxdata_grp.name.endswith("hplc")  # we are reading HPLC data
            signal = nxdata_grp.attrs["signal"]
            axis = nxdata_grp.attrs["axes"]
            Isum = nxdata_grp[signal][()]
            idx = nxdata_grp[axis][()]
            integrated = nxdata_grp.parent["results"]
            signal = integrated.attrs["signal"]
            I = integrated[signal][()]
            axes = integrated.attrs["axes"][-1]
            q = integrated[axes][()]
            sigma = integrated["errors"][()]

            npt = len(q)
            unit = pyFAI.units.to_unit(axes + "_" + integrated[axes].attrs["units"])
            integration_grp = nxdata_grp.parent
            poni = str(integration_grp["configuration/file_name"][()]).strip()
            if not os.path.exists(poni):
                poni = str(integration_grp["configuration/data"][()]).strip()
            polarization = integration_grp["configuration/polarization_factor"][()]
            method = IntegrationMethod.select_method(**json.loads(integration_grp["configuration/integration_method"][()]))[0]
            instrument_grp = nxsr.get_class(entry_grp, class_type="NXinstrument")[0]
            detector_grp = nxsr.get_class(instrument_grp, class_type="NXdetector")[0]
            mask = detector_grp["pixel_mask"].attrs["filename"]
            mono_grp = nxsr.get_class(instrument_grp, class_type="NXmonochromator")[0]
            energy = mono_grp["energy"][()]

            # Read the sample description:
            sample_grp = nxsr.get_class(entry_grp, class_type="NXsample")[0]
            sample_name = posixpath.split(sample_grp.name)[-1]

            buffer = sample_grp["buffer"][()] if "buffer" in sample_grp else ""
            concentration = sample_grp["concentration"][()] if "concentration" in sample_grp else ""
            description = sample_grp["description"][()] if "description" in sample_grp else ""
            hplc = sample_grp["hplc"][()] if "hplc" in sample_grp else ""
            temperature = sample_grp["temperature"][()] if "temperature" in sample_grp else ""
            temperature_env = sample_grp["temperature_env"][()] if "temperature_env" in sample_grp else ""
            sample = Sample(sample_name, description, buffer, concentration, hplc, temperature_env, temperature)
            meas_grp = nxsr.get_class(entry_grp, class_type="NXdata")[0]
            timestamps = []
            for ts_name in ("timestamps", "time-stamps"):
                if ts_name in meas_grp:
                    timestamps = meas_grp[ts_name][()]
                    break
        return NexusJuice(filename, h5path, npt, unit, idx, Isum, q, I, sigma, poni, mask, energy, polarization, method, sample, timestamps)
        "filename h5path npt unit idx Isum q I sigma poni mask energy polarization method sample timestamps"
