#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* HPLC mode: Rebuild the complete chromatogram and perform basic analysis on it.
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "14/04/2021"
__status__ = "development"
__version__ = "0.1.0"

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
import scipy.signal
import scipy.ndimage
from sklearn.decomposition import NMF
from .common import Sample, Ispyb, get_equivalent_frames, cmp_float, get_integrator, KeyCache, \
                    polarization_factor, method, Nexus, get_isotime, SAXS_STYLE, NORMAL_STYLE, \
                    Sample, create_nexus_sample
from .ispyb import IspybConnector

NexusJuice = namedtuple("NexusJuice", "filename h5path npt unit idx Isum q I sigma poni mask energy polarization method sample")


def smooth_chromatogram(signal, window):
    """smooth-out the chromatogram
    
    :param signal: the chomatogram as 1d array
    :param window: the size of the window
    """
    # wodd = window+1 if window%2==0 else window
    wodd = 2 * window + 1
    smth = scipy.signal.medfilt(signal, wodd)
    # sigma = wmin/2*numpy.sqrt(2*numpy.log(2))
    # w2 = int(6*sigma)
    # if w2%2 == 0: w2+=1
    # print(wmin, sigma, w2)
#     g = scipy.signal.gaussian(wodd, wodd/4)
#     g /= g.sum()

    # small kernel smoothing to remove steps induced by medfilt
    # smth2 = scipy.signal.convolve(signal, g,"same")
    return smth


def search_peaks(signal, wmin=10, scale=0.9):
    """
    Label all peak regions of chromatogram.
    
    :param signal=smooth signal
    :param wmin: minimum width for a peak. smaller ones are discarded.
    :param scale: shrink factor (i.e. <1 for the search zone)
    """

    smth = smooth_chromatogram(signal, window=wmin)
    res = numpy.zeros(signal.shape, dtype=numpy.uint8)
    w = signal.size
    while w > wmin:
        peaks = scipy.signal.find_peaks_cwt(smth, [w])
        if len(peaks):
            widths = scipy.signal.peak_widths(smth, peaks)
            m = widths[0] >= wmin
            if m.max():
                for p, q in zip(peaks[m], widths[0][m]):
                    # print(w, p, q)
                    q = int(numpy.ceil(q))
                    if q > p:
                        start = smth[:q]
                    else:
                        start = smth[p - q:p]
                    if p + q >= signal.size:
                        stop = smth[-q:]
                    else:
                        stop = smth[p:q + p]
                    pos = numpy.argmin(abs(start - stop))
                    res[p + pos - q: p + pos] = 1
        w *= scale
    return scipy.ndimage.label(res)


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
       "nmf_components = 5, 
      "wait_for": [jobid_img001, jobid_img002],
      "plugin_name": "bm29.hplc"
    } 
    """
    NMF_COMP = 5
    "Default number of Non-negative matrix factorisation components. Correspond to the number of spieces"

    def __init__(self):
        Plugin.__init__(self)
        self.input_files = []
        self.sample_file = None
        self.nxs = None
        self.output_file = None
        self.juices = []
        self.nmf_components = self.NMF_COMP
        self.to_pyarch = {}
        self.ispyb = None

    def setup(self):
        Plugin.setup(self)

        for job_id in self.input.get("wait_for", []):
            self.wait_for(job_id)

        self.input_files = [os.path.abspath(i) for i in self.input.get("integrated_files", "")]

        self.output_file = self.input.get("output_file")
        if not self.output_file:
            self.output_file = os.path.commonprefix(self.input_files) + "_hplc.h5"
            self.log_warning("No output file provided, using " + self.output_file)
        self.nmf_components = int(self.input.get("nmf_components", self.NMF_COMP))
        self.ispyb = Ispyb._fromdict(self.input.get("ispyb", {}))

    def process(self):
        self.create_nexus()
        if not self.input.get("no_ispyb"):
            self.send_to_ispyb()

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

        q = self.juices[0].q
        unit = self.juices[0].unit
        radial_unit, unit_name = str(unit).split("_", 1)

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

    # Process 3: NMF matrix decomposition
        nmf_grp = nxs.new_class(entry_grp, "3_NMF", "NXprocess")
        nmf_grp["sequence_index"] = 3
        nmf = NMF(n_components=self.nmf_components, init='nndsvd',
                  max_iter=1000)
        W = nmf.fit_transform(I.T)
        eigen_data = nxs.new_class(nmf_grp, "eigenvectors", "NXdata")
        eigen_ds = eigen_data.create_dataset("W", data=numpy.ascontiguousarray(W.T, dtype=numpy.float32))
        eigen_ds.attrs["interpretation"] = "spectrum"
        eigen_data.attrs["signal"] = "W"
        eigen_data.attrs["SILX_style"] = SAXS_STYLE

        eigen_ds.attrs["units"] = "arbitrary"
        eigen_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"

        H = nmf.components_
        chroma_data = nxs.new_class(nmf_grp, "chromatogram", "NXdata")
        chroma_ds = chroma_data.create_dataset("H", data=numpy.ascontiguousarray(H, dtype=numpy.float32))
        chroma_ds.attrs["interpretation"] = "spectrum"
        chroma_data.attrs["signal"] = "H"

        quantiles = (0.1, 0.6)  # assume weakest diffracting are background keep 10-40%
        background = W[:, 0] * (numpy.sort(H[0])[int(nframes * quantiles[0]): int(nframes * quantiles[-1])]).mean()
        bg_data = nxs.new_class(nmf_grp, "background", "NXdata")
        bg_data.attrs["signal"] = "I"
        bg_data.attrs["SILX_style"] = SAXS_STYLE
        bg_data.attrs["axes"] = radial_unit
        bg_ds = bg_data.create_dataset("I", data=numpy.ascontiguousarray(background, dtype=numpy.float32))
        bg_ds.attrs["interpretation"] = "spectrum"
        bg_data.attrs["quantiles"] = quantiles
        bg_q_ds = bg_data.create_dataset(radial_unit,
                                         data=numpy.ascontiguousarray(q, dtype=numpy.float32))
        bg_q_ds.attrs["units"] = unit_name
        radius_unit = "nm" if "nm" in unit_name else "Å"
        bg_q_ds.attrs["long_name"] = f"Scattering vector q ({radius_unit}⁻¹)"

    # Process 4: fraction of chromatogram analysis
        fraction_grp = nxs.new_class(entry_grp, "4_fractions", "NXprocess")
        fraction_grp["sequence_index"] = 4
        fraction_grp["width"] = window = 10

        I_sub = I - background
        fractions, nfractions = search_peaks(Isum, window)
        if nfractions:
            for fraction in scipy.ndimage.find_objects(fractions, nfractions):
                self.one_fraction(fraction[0], nxs, fraction_grp, I_sub, sigma)

    def one_fraction(self, fraction, nxs, top_grp, I_sub, sigma):
        """
        :param fraction: slice with start and end
        :param nxs: opened Nexus file object
        
        """
        q = self.juices[0].q
        unit = self.juices[0].unit
        sample = self.juices[0].sample
        radial_unit, unit_name = str(unit).split("_", 1)

        f_grp = nxs.new_class(top_grp, f"{fraction.start}-{fraction.stop}", "NXprocess")
        f_grp["sequence_index"] = 5
        f_grp["program"] = "dahu.plugins.bm29.hplc"
        f_grp["version"] = __version__
        f_grp["date"] = get_isotime()
        
        avg_data = nxs.new_class(f_grp, "1_average", "NXdata")
        avg_data.attrs["SILX_style"] = SAXS_STYLE
        avg_data.attrs["title"] = f"{sample.name}, frames {fraction.start}-{fraction.stop} averaged, buffer subtracted"
        avg_data.attrs["signal"] = "I"
        avg_data.attrs["axes"] = radial_unit
        f_grp.attrs["default"] = avg_data.name
        avg_q_ds = avg_data.create_dataset(radial_unit,
                                           data=numpy.ascontiguousarray(q, dtype=numpy.float32))
        avg_q_ds.attrs["units"] = unit_name
        radius_unit = "nm" if "nm" in unit_name else "Å"
        avg_q_ds.attrs["long_name"] = f"Scattering vector q ({radius_unit}⁻¹)"
        I_frc = I_sub[fraction].mean(axis=0)
        fsig2 = sigma[fraction] ** 2
        sigma_frc = numpy.sqrt(fsig2.sum(axis=0) / fsig2.shape[0])
        ai2_int_ds = avg_data.create_dataset("I", data=numpy.ascontiguousarray(I_frc, dtype=numpy.float32))
        ai2_std_ds = avg_data.create_dataset("errors",
                                             data=numpy.ascontiguousarray(sigma_frc, dtype=numpy.float32))

        ai2_int_ds.attrs["interpretation"] = "spectrum"
        ai2_int_ds.attrs["units"] = "arbitrary"
        ai2_int_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"
        #  ai2_int_ds.attrs["uncertainties"] = "errors" #this does not work
        ai2_std_ds.attrs["interpretation"] = "spectrum"
        ai2_int_ds.attrs["units"] = "arbitrary"


    # Process 4: Guinier analysis
        guinier_grp = nxs.new_class(f_grp, "6_Guinier_analysis", "NXprocess")
        guinier_grp["sequence_index"] = 6
        guinier_grp["program"] = "freesas.autorg"
        guinier_grp["version"] = freesas.version
        guinier_grp["date"] = get_isotime()
        guinier_autorg = nxs.new_class(guinier_grp, "autorg", "NXcollection")
        guinier_gpa = nxs.new_class(guinier_grp, "gpa", "NXcollection")
        guinier_guinier = nxs.new_class(guinier_grp, "guinier", "NXcollection")
        guinier_data = nxs.new_class(guinier_grp, "results", "NXdata")
        guinier_data.attrs["SILX_style"] = NORMAL_STYLE
        guinier_data.attrs["title"] = "Guinier analysis"
    # Stage4 processing: autorg and auto_gpa
        sasm = numpy.vstack((q, I_frc, sigma_frc)).T

        try:
            gpa = auto_gpa(sasm)
        except Exception as error:
            guinier_gpa["Failed"] = "%s: %s" % (error.__class__.__name__, error)
            gpa = None
        else:
            #  "Rg sigma_Rg I0 sigma_I0 start_point end_point quality aggregated"
            guinier_gpa["Rg"] = gpa.Rg
            guinier_gpa["Rg"].attrs["unit"] = radius_unit
            guinier_gpa["Rg_error"] = gpa.sigma_Rg
            guinier_gpa["Rg_error"].attrs["unit"] = radius_unit
            guinier_gpa["I0"] = gpa.I0
            guinier_gpa["I0_error"] = gpa.sigma_I0
            guinier_gpa["start_point"] = gpa.start_point
            guinier_gpa["end_point"] = gpa.end_point
            guinier_gpa["quality"] = gpa.quality
            guinier_gpa["aggregated"] = gpa.aggregated

        try:
            guinier = auto_guinier(sasm)
        except Exception as error:
            guinier_guinier["Failed"] = "%s: %s" % (error.__class__.__name__, error)
            guinier = None
        else:
            #  "Rg sigma_Rg I0 sigma_I0 start_point end_point quality aggregated"
            guinier_guinier["Rg"] = guinier.Rg
            guinier_guinier["Rg"].attrs["unit"] = radius_unit
            guinier_guinier["Rg_error"] = guinier.sigma_Rg
            guinier_guinier["Rg_error"].attrs["unit"] = radius_unit
            guinier_guinier["I0"] = guinier.I0
            guinier_guinier["I0_error"] = guinier.sigma_I0
            guinier_guinier["start_point"] = guinier.start_point
            guinier_guinier["end_point"] = guinier.end_point
            guinier_guinier["quality"] = guinier.quality
            guinier_guinier["aggregated"] = guinier.aggregated
            guinier_guinier["qₘᵢₙ·Rg"] = guinier.Rg * q[guinier.start_point]
            guinier_guinier["qₘₐₓ·Rg"] = guinier.Rg * q[guinier.end_point - 1]

        try:
            autorg = autoRg(sasm)
        except Exception as err:
            guinier_autorg["Failed"] = "%s: %s" % (err.__class__.__name__, err)
            autorg = None
        else:
            if autorg.Rg < 0:
                guinier_autorg["Failed"] = "No Guinier region found with this algorithm"
                autorg = None
            else:
                # "Rg sigma_Rg I0 sigma_I0 start_point end_point quality aggregated"
                guinier_autorg["Rg"] = autorg.Rg
                guinier_autorg["Rg"].attrs["unit"] = radius_unit
                guinier_autorg["Rg_error"] = autorg.sigma_Rg
                guinier_autorg["Rg_error"].attrs["unit"] = radius_unit
                guinier_autorg["I0"] = autorg.I0
                guinier_autorg["I0_error"] = autorg.sigma_I0
                guinier_autorg["start_point"] = autorg.start_point
                guinier_autorg["end_point"] = autorg.end_point
                guinier_autorg["quality"] = autorg.quality
                guinier_autorg["aggregated"] = autorg.aggregated
                guinier_autorg["qₘᵢₙ·Rg"] = autorg.Rg * q[autorg.start_point]
                guinier_autorg["qₘₐₓ·Rg"] = autorg.Rg * q[autorg.end_point - 1]

        #  take one of the fits
        if guinier:
            guinier_data["source"] = "auto_guinier"
        elif autorg:
            guinier = autorg
            guinier_data["source"] = "autorg"
        elif gpa:
            guinier = gpa
            guinier_data["source"] = "gpa"
        else:
            guinier = None
            guinier_data["source"] = "None"

    # Stage #4 Guinier plot generation:

        q, I, err = sasm.T[:3]
        mask = (I > 0) & numpy.isfinite(I) & (q > 0) & numpy.isfinite(q)
        if err is not None:
            mask &= (err > 0.0) & numpy.isfinite(err)
        mask = mask.astype(bool)
        if guinier:
            intercept = numpy.log(guinier.I0)
            slope = -guinier.Rg ** 2 / 3.0
            end = numpy.where(q > 1.5 / guinier.Rg)[0][0]
            mask[end:] = False

        q2 = q[mask] ** 2
        logI = numpy.log(I[mask])
        dlogI = err[mask] / logI
        q2_ds = guinier_data.create_dataset("q2", data=q2.astype(numpy.float32))
        q2_ds.attrs["unit"] = radius_unit + "⁻²"
        q2_ds.attrs["long_name"] = "q² (%s⁻²)" % radius_unit
        q2_ds.attrs["interpretation"] = "spectrum"
        lnI_ds = guinier_data.create_dataset("logI", data=logI.astype(numpy.float32))
        lnI_ds.attrs["long_name"] = "log(I)"
        lnI_ds.attrs["interpretation"] = "spectrum"
        erI_ds = guinier_data.create_dataset("errors", data=dlogI.astype(numpy.float32))
        erI_ds.attrs["interpretation"] = "spectrum"

        if guinier:
            guinier_data["fit"] = intercept + slope * q2
            guinier_data["fit"].attrs["slope"] = slope
            guinier_data["fit"].attrs["intercept"] = intercept

        guinier_data_attrs = guinier_data.attrs
        guinier_data_attrs["signal"] = "logI"
        guinier_data_attrs["axes"] = "q2"
        guinier_data_attrs["auxiliary_signals"] = "fit"
        guinier_grp.attrs["default"] = guinier_data.name
        if guinier is None:
            f_grp.attrs["default"] = avg_data.name
            self.log_error("No Guinier region found, data of dubious quality", do_raise=False)
            return

    # Process 5: Kratky plot
        kratky_grp = nxs.new_class(f_grp, "7_dimensionless_Kratky_plot", "NXprocess")
        kratky_grp["sequence_index"] = 7
        kratky_grp["program"] = "freesas.autorg"
        kratky_grp["version"] = freesas.version
        kratky_grp["date"] = get_isotime()
        kratky_data = nxs.new_class(kratky_grp, "results", "NXdata")
        kratky_data.attrs["SILX_style"] = NORMAL_STYLE
        kratky_data.attrs["title"] = "Dimensionless Kratky plots"
        kratky_grp.attrs["default"] = kratky_data.name

    # Stage #5 Kratky plot generation:
        Rg = guinier.Rg
        I0 = guinier.I0
        xdata = q * Rg
        ydata = xdata * xdata * I / I0
        dy = xdata * xdata * err / I0
        qRg_ds = kratky_data.create_dataset("qRg", data=xdata.astype(numpy.float32))
        qRg_ds.attrs["interpretation"] = "spectrum"
        qRg_ds.attrs["long_name"] = "q·Rg (unit-less)"
        k_ds = kratky_data.create_dataset("q2Rg2I/I0", data=ydata.astype(numpy.float32))
        k_ds.attrs["interpretation"] = "spectrum"
        k_ds.attrs["long_name"] = "q²Rg²I(q)/I₀"
        ke_ds = kratky_data.create_dataset("errors", data=dy.astype(numpy.float32))
        ke_ds.attrs["interpretation"] = "spectrum"
        kratky_data_attrs = kratky_data.attrs
        kratky_data_attrs["signal"] = "q2Rg2I/I0"
        kratky_data_attrs["axes"] = "qRg"

    # stage 6: Rambo-Tainer invariant
        rti_grp = nxs.new_class(f_grp, "8_invariants", "NXprocess")
        rti_grp["sequence_index"] = 8
        rti_grp["program"] = "freesas.invariants"
        rti_grp["version"] = freesas.version
        rti_data = nxs.new_class(rti_grp, "results", "NXdata")
        # average_data.attrs["SILX_style"] = SAXS_STYLE
        # average_data.attrs["signal"] = "intensity_normed"
        # Rambo_Tainer
        rti = freesas.invariants.calc_Rambo_Tainer(sasm, guinier)
        Vc_ds = rti_data.create_dataset("Vc", data=rti.Vc)
        Vc_ds.attrs["unit"] = "nm²"
        Vc_ds.attrs["formula"] = "Rambo-Tainer: Vc = I₀/(sum_q qI(q) dq)"
        sigma_Vc_ds = rti_data.create_dataset("Vc_error", data=rti.sigma_Vc)
        sigma_Vc_ds.attrs["unit"] = "nm²"

        Qr_ds = rti_data.create_dataset("Qr", data=rti.Qr)
        Qr_ds.attrs["unit"] = "nm"
        Qr_ds.attrs["formula"] = "Rambo-Tainer: Qr = Vc/Rg"
        sigma_Qr_ds = rti_data.create_dataset("Qr_error", data=rti.sigma_Qr)
        sigma_Qr_ds.attrs["unit"] = "nm"

        mass_ds = rti_data.create_dataset("mass", data=rti.mass)
        mass_ds.attrs["unit"] = "kDa"
        mass_ds.attrs["formula"] = "Rambo-Tainer: mass = (Qr/ec)^(1/k)"
        sigma_mass_ds = rti_data.create_dataset("mass_error", data=rti.sigma_mass)
        sigma_mass_ds.attrs["unit"] = "kDa"

        self.Vc = rti.Vc
        self.mass = rti.mass

        volume = self.to_pyarch["volume"] = freesas.invariants.calc_Porod(sasm, guinier)
        volume_ds = rti_data.create_dataset("volume", data=volume)
        volume_ds.attrs["unit"] = "nm³"
        volume_ds.attrs["formula"] = "Porod: V = 2*π²I₀²/(sum_q I(q)q² dq)"
        self.to_pyarch["rti"] = rti

    # stage 7: Pair distribution function, what is the equivalent of datgnom
        bift_grp = nxs.new_class(f_grp, "9_indirect_Fourier_transformation", "NXprocess")
        bift_grp["sequence_index"] = 9
        bift_grp["program"] = "freesas.bift"
        bift_grp["version"] = freesas.version
        bift_grp["date"] = get_isotime()
        bift_data = nxs.new_class(bift_grp, "results", "NXdata")
        bift_data.attrs["SILX_style"] = NORMAL_STYLE
        bift_data.attrs["title"] = "Pair distance distribution function p(r)"

        cfg_grp = nxs.new_class(bift_grp, "configuration", "NXcollection")
    # Process stage7, i.e. perform the IFT
        try:
            bo = BIFT(q, I, err)
            cfg_grp["Rg"] = guinier.Rg
            # Pretty limited quality as we have real time constrains
            cfg_grp["npt"] = npt = 64
            cfg_grp["Dmax÷Rg"] = 3
            Dmax = bo.set_Guinier(guinier, Dmax_over_Rg=3)
            # Pretty limited quality as we have real time constrains

            # First scan on alpha:
            cfg_grp["alpha_sup"] = alpha_max = bo.guess_alpha_max(npt)
            cfg_grp["alpha_inf"] = 1 / alpha_max
            cfg_grp["alpha_scan_steps"] = 11

            key = bo.grid_scan(Dmax, Dmax, 1,
                               1.0 / alpha_max, alpha_max, 11, npt)
            Dmax, alpha = key[:2]
            # Then scan on Dmax:
            cfg_grp["Dmax_sup"] = guinier.Rg * 4
            cfg_grp["Dmax_inf"] = guinier.Rg * 2
            cfg_grp["Dmax_scan_steps"] = 5
            key = bo.grid_scan(guinier.Rg * 2, guinier.Rg * 4, 5,
                               alpha, alpha, 1, npt)
            Dmax, alpha = key[:2]
            if bo.evidence_cache[key].converged:
                bo.update_wisdom()
                use_wisdom = True
            else:
                use_wisdom = False
            res = minimize(bo.opti_evidence, (Dmax, log(alpha)), args=(npt, use_wisdom), method="powell")
            cfg_grp["Powell_steps"] = res.nfev
            cfg_grp["Monte-Carlo_steps"] = 0
        except Exception as error:
            bift_grp["Failed"] = "%s: %s" % (error.__class__.__name__, error)
            bo = None
        else:
            stats = bo.calc_stats()
            bift_grp["alpha"] = stats.alpha_avg
            bift_grp["alpha_error"] = stats.alpha_std
            self.Dmax = bift_grp["Dmax"] = stats.Dmax_avg
            bift_grp["Dmax_error"] = stats.Dmax_std
            bift_grp["S0"] = stats.regularization_avg
            bift_grp["S0_error"] = stats.regularization_std
            bift_grp["Chi2r"] = stats.chi2r_avg
            bift_grp["Chi2r_error"] = stats.chi2r_std
            bift_grp["logP"] = stats.evidence_avg
            bift_grp["logP_error"] = stats.evidence_std
            bift_grp["Rg"] = stats.Rg_avg
            bift_grp["Rg_error"] = stats.Rg_std
            bift_grp["I0"] = stats.I0_avg
            bift_grp["I0_error"] = stats.I0_std
            # Now the plot:
            r_ds = bift_data.create_dataset("r", data=stats.radius.astype(numpy.float32))
            r_ds.attrs["interpretation"] = "spectrum"

            r_ds.attrs["unit"] = radius_unit
            r_ds.attrs["long_name"] = "radius r(%s)" % radius_unit
            p_ds = bift_data.create_dataset("p(r)", data=stats.density_avg.astype(numpy.float32))
            p_ds.attrs["interpretation"] = "spectrum"
            bift_data["errors"] = stats.density_std
            bift_data.attrs["signal"] = "p(r)"
            bift_data.attrs["axes"] = "r"

            r = stats.radius
            T = numpy.outer(q, r / pi)
            T = (4 * pi * (r[-1] - r[0]) / (len(r) - 1)) * numpy.sinc(T)
            bift_ds = avg_data.create_dataset("BIFT", data=T.dot(stats.density_avg).astype(numpy.float32))
            bift_ds.attrs["interpretation"] = "spectrum"
            avg_data.attrs["auxiliary_signals"] = "BIFT"
            bift_grp.attrs["default"] = bift_data.name

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

    def send_to_ispyb(self):
        """Data sent to ISPyB are:
            * hdf5File
            * jsonFile built from HDF5
            * hplcPlot various plots generated
        """
        if self.ispyb and self.ispyb.url and parse_url(self.ispyb.url).host:
            ispyb = IspybConnector(*self.ispyb)

            # Composition of to_pyarch:
            # X self.to_pyarch["id"] = self.id
            # self.to_pyarch["buffer"] = self.buffer
            # self.to_pyarch["first_curve"] = self.first_curve
            # self.to_pyarch["frames"] = {}
            # self.to_pyarch["curves"] = self.curves
            # self.to_pyarch["for_buffer"] = self.for_buffer
            # X self.to_pyarch["hdf5_filename"] = self.hdf5_filename
            # X self.to_pyarch["chunk_size"] = self.chunk_size

            self.to_pyarch["hdf5_filename"] = self.output_file
            self.to_pyarch["chunk_size"] = min(j.Isum.size for j in self.some_juice)
            self.to_pyarch["id"] = os.path.commonprefix(self.input_files)
            ispyb.send_hplc(self.to_pyarch)
        else:
            self.log_warning(f"Not sending to ISPyB: no valid URL in {self.ispyb}")
