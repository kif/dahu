#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* HPLC mode: Rebuild the complete chromatogram and perform basic analysis on it.
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "07/10/2020"
__status__ = "development"
__version__ = "0.0.1"

import os
import json
from math import log, pi
from collections import namedtuple
from urllib3.util import parse_url
from dahu.plugin import Plugin
from dahu.utils import fully_qualified_name
import logging
logger = logging.getLogger("bm29.hplc")
import numpy
import h5py
import pyFAI, pyFAI.azimuthalIntegrator
from pyFAI.method_registry import IntegrationMethod
import freesas, freesas.cormap, freesas.invariants
from freesas.autorg import auto_gpa, autoRg, auto_guinier
from freesas.bift import BIFT
from scipy.optimize import minimize
from .common import Sample, Ispyb, get_equivalent_frames, cmp_float, get_integrator, KeyCache, \
                    polarization_factor, method, Nexus, get_isotime, SAXS_STYLE, NORMAL_STYLE, \
                    Sample, create_nexus_sample
from .ispyb import IspybConnector

NexusJuice = namedtuple("NexusJuice", "filename h5path npt unit idx Isum q I sigma poni mask energy polarization method sample")


class HPLC(Plugin):
    """ Rebuild the complete chromatogram and perform basic analysis on it.
    
        Typical JSON file:
    {
      "integrated_files": ["img_001.h5", "img_002.h5"],
      "output_file": "hplc.h5"
      "ispyb": {
        "url": "http://ispyb.esrf.fr:1234",
        "pyarch": "/data/pyarch/mx1234/sample", 
        "measurement_id": -1,
        "collection_id": -1
       },
      "wait_for": [jobid_img001, jobid_img002],
      "plugin_name": "bm29.hplc"
    } 
    """

    def __init__(self):
        Plugin.__init__(self)
        self.input_files = []
        self.sample_file = None
        self.nxs = None
        self.output_file = None
        self.juices = []

    def setup(self):
        Plugin.setup(self)
        self.input_files = [os.path.abspath(i) for i in self.input.get("integrated_files", "")]

        self.output_file = self.input.get("output_file")
        if self.output_file is None:
            self.output_file = os.path.commonprefix(self.input_files) + "_hplc.h5"
            self.log_warning("No output file provided, using " + self.output_file)

    def process(self):
        self.some_juice = [self.read_nexus(i) for i in self.input_files]
        self.create_nexus()

    def teardown(self):
        Plugin.teardown(self)

    def create_nexus(self):
        nxs = self.nxs = Nexus(self.output_file, mode="w")
        entry_grp = nxs.new_entry("entry", self.input.get("plugin_name", "dahu"),
                              title='BioSaxs HPLC experiment',
                              force_time=get_isotime())
        nxs.h5.attrs["default"] = entry_grp.name

    # Configuration
        cfg_grp = nxs.new_class(entry_grp, "configuration", "NXnote")
        cfg_grp.create_dataset("data", data=json.dumps(self.input, indent=2, separators=(",\r\n", ":\t")))
        cfg_grp.create_dataset("format", data="text/json")

    # Process 0: Measurement group
        input_grp = nxs.new_class(entry_grp, "0_measurement", "NXcollection")
        input_grp["sequence_index"] = 0

        for idx, filename in enumerate(self.input_files):
            juice = self.read_nexus(filename)
            if juice is not None:
                rel_path = os.path.relpath(os.path.abspath(filename), os.path.dirname(os.path.abspath(self.output_file)))
                input_grp["LImA_%04i" % idx] = h5py.ExternalLink(rel_path, juice.h5path)
                self.juices.append(juice)

        # Sample: outsourced !
        create_nexus_sample(nxs, entry_grp, self.juices[0].sample)

    # Process 1: Chromatogram
        chroma_grp = nxs.new_class(entry_grp, "1_chromatogram", "NXprocess")
        chroma_grp["sequence_index"] = 1
        nframes = max(i.idx.max() for i in self.juices) + 1
        nbin = self.juices[0].q.size

        I = numpy.zeros((nframes, nbin), dtype=numpy.float32)
        sigma = numpy.zeros((nframes, nbin), dtype=numpy.float32)
        Isum = numpy.zeros(nframes)

        ids = numpy.arange(nframes)
        idx = numpy.concatenate([i.idx for i in self.juices])
        I[idx] = numpy.vstack([i.I for i in self.juices])
        Isum[idx] = numpy.concatenate([i.Isum for i in self.juices])
        sigma[idx] = numpy.vstack([i.sigma for i in self.juices])

        hplc_data = nxs.new_class(chroma_grp, "hplc", "NXdata")
        hplc_data.attrs["title"] = "Chromatogram"
        sum_ds = hplc_data.create_dataset("sum", data=Isum, dtype=numpy.float32)
        sum_ds.attrs["interpretation"] = "spectrum"
        sum_ds.attrs["long_name"] = "Summed Intensity"
        frame_ds = hplc_data.create_dataset("frame_ids", data=ids, dtype=numpy.uint32)
        frame_ds.attrs["interpretation"] = "spectrum"
        frame_ds.attrs["long_name"] = "frame index"
        hplc_data.attrs["signal"] = "sum"
        hplc_data.attrs["axes"] = "frame_ids"
        entry_grp.attrs["default"] = hplc_data.name
        integration_data = nxs.new_class(chroma_grp, "results", "NXdata")
        chroma_grp.attrs["title"] = str(self.juices[0].sample)

        int_ds = integration_data.create_dataset("I", data=numpy.ascontiguousarray(I, dtype=numpy.float32))
        std_ds = integration_data.create_dataset("errors", data=numpy.ascontiguousarray(sigma, dtype=numpy.float32))
        q_ds = integration_data.create_dataset("q", data=self.juices[0].q)
        q_ds.attrs["interpretation"] = "spectrum"
        integration_data.attrs["signal"] = "I"
        integration_data.attrs["axes"] = [".", "q"]
        integration_data.attrs["SILX_style"] = SAXS_STYLE

        int_ds.attrs["interpretation"] = "spectrum"
        int_ds.attrs["units"] = "arbitrary"
        int_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"
        # int_ds.attrs["uncertainties"] = "errors" This does not work
        int_ds.attrs["scale"] = "log"
        std_ds.attrs["interpretation"] = "spectrum"

    # Process 2: SVD decomposition
        svd_grp = nxs.new_class(entry_grp, "2_SVD", "NXprocess")
        svd_grp["sequence_index"] = 2
        logi = numpy.arcsinh(I.T)
        U, S, V = numpy.linalg.svd(logi, full_matrices=False)

        # Number of Eignevector to keep:
        svd_grp["Ref"] = "https://arxiv.org/pdf/1305.5870.pdf"
        beta = nframes / nbin if nframes <= nbin else 1.0
        omega = 0.56 * beta ** 3 - 0.95 * beta ** 2 + 1.82 * beta + 1.43
        tau = numpy.median(S) * omega
        r = numpy.sum(S > tau)

        # Flip axis with negative signal
        flip = V.max(axis=1) < -V.min(axis=1)
        nflip = numpy.where(flip)
        V[nflip] = -V[nflip]
        U[:, nflip] = -U[:, nflip]

        eigen_data = nxs.new_class(svd_grp, "eigenvectors", "NXdata")
        eigen_ds = eigen_data.create_dataset("U", data=numpy.ascontiguousarray(U.T[:r], dtype=numpy.float32))
        eigen_ds.attrs["interpretation"] = "spectrum"
        eigen_data.attrs["signal"] = "U"

        chroma_data = nxs.new_class(svd_grp, "chromatogram", "NXdata")
        chroma_ds = chroma_data.create_dataset("V", data=numpy.ascontiguousarray(V[:r], dtype=numpy.float32))
        chroma_ds.attrs["interpretation"] = "spectrum"
        chroma_data.attrs["signal"] = "V"

        svd_grp.create_dataset("eigenvalues", data=S[:r], dtype=numpy.float32)

    @staticmethod
    def read_nexus(filename):
        "return some NexusJuice from a HDF5 file "
        with Nexus(filename, "r") as nxsr:
            entry_name = nxsr.h5.attrs["default"]
            entry_grp = nxsr.h5[entry_name]
            h5path = entry_grp.name
            nxdata_grp = nxsr.h5[entry_grp.attrs["default"]]
            assert nxdata_grp.name.endswith("hplc")  # we are reading HPLC data
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
#             img_grp = nxsr.get_class(entry_grp["3_time_average"], class_type="NXdata")[0]
#             image2d = img_grp["intensity_normed"][()]
#             error2d = img_grp["intensity_std"][()]
            # Read the sample description:
            sample_grp = nxsr.get_class(entry_grp, class_type="NXsample")[0]
            sample_name = sample_grp.name

            buffer = sample_grp["buffer"][()] if "buffer" in sample_grp else ""
            concentration = sample_grp["concentration"][()] if "concentration" in sample_grp else ""
            description = sample_grp["description"][()] if "description" in sample_grp else ""
            hplc = sample_grp["hplc"][()] if "hplc" in sample_grp else ""
            temperature = sample_grp["temperature"][()] if "temperature" in sample_grp else ""
            temperature_env = sample_grp["temperature_env"][()] if "temperature_env" in sample_grp else ""
            sample = Sample(sample_name, description, buffer, concentration, hplc, temperature_env, temperature)

        return NexusJuice(filename, h5path, npt, unit, idx, Isum, q, I, sigma, poni, mask, energy, polarization, method, sample)
        "filename h5path npt unit idx Isum q I sigma poni mask energy polarization method sample"
