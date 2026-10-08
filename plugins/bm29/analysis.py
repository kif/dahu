"""Data Analysis for BM29: BioSaxs

Common SAXS analysis sequence shared by the SubtractBuffer and the HPLC plugins:

* Guinier analysis (auto_guinier, autorg and GPA)
* Dimensionless Kratky plot
* Invariants (Rambo-Tainer and Porod)
* Indirect Fourier transformation (BIFT)

Each step stores its results as a NXprocess in a Nexus file.
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "07/10/2026"
__status__ = "development"
__version__ = "0.1.0"

import itertools
import logging
import posixpath
from math import log, pi
from typing import NamedTuple

import freesas
import freesas.invariants
import numpy
from freesas.autorg import auto_gpa, auto_guinier, autoRg
from freesas.bift import BIFT
from scipy.optimize import minimize

from .common import NORMAL_STYLE
from .nexus import get_isotime

logger = logging.getLogger("bm29.analysis")


class AnalysisResult(NamedTuple):
    """Results of the SAXS analysis sequence. Fields are None when the step failed or was not reached"""
    guinier: object = None  # RG_RESULT of the selected Guinier fit
    rti: object = None  # Rambo-Tainer invariants
    volume: float = None  # Porod volume
    bift: object = None  # BIFT statistics (StatsResult)


def _new_process(nxs, parent_grp, name, sequence_index, program):
    "Create a NXprocess group with the common metadata"
    grp = nxs.new_class(parent_grp, name, "NXprocess")
    grp["sequence_index"] = sequence_index
    grp["program"] = program
    grp["version"] = freesas.version
    grp["date"] = get_isotime()
    return grp


def _store_rg(grp, fit, radius_unit, q=None):
    "Store a RG_RESULT into a NXcollection, with qmin·Rg and qmax·Rg if q is provided"
    #  "Rg sigma_Rg I0 sigma_I0 start_point end_point quality aggregated"
    grp["Rg"] = fit.Rg
    grp["Rg"].attrs["unit"] = radius_unit
    grp["Rg_error"] = fit.sigma_Rg
    grp["Rg_error"].attrs["unit"] = radius_unit
    grp["I0"] = fit.I0
    grp["I0_error"] = fit.sigma_I0
    grp["start_point"] = fit.start_point
    grp["end_point"] = fit.end_point
    grp["quality"] = fit.quality
    grp["aggregated"] = fit.aggregated
    if q is not None:
        grp["qₘᵢₙ·Rg"] = fit.Rg * q[fit.start_point]
        grp["qₘₐₓ·Rg"] = fit.Rg * q[fit.end_point - 1]


def guinier_analysis(nxs, parent_grp, name, sequence_index, sasm, radius_unit):
    """Perform the Guinier analysis with 3 algorithms (GPA, auto_guinier, autorg) and keep the best one.

    :param nxs: opened Nexus file object
    :param parent_grp: group in which the NXprocess is created
    :param name: name of the NXprocess group
    :param sequence_index: index of the process
    :param sasm: 2D array with q, I, sigma as columns
    :param radius_unit: "nm" or "Å"
    :return: the selected Guinier fit (RG_RESULT) or None if all algorithms failed
    """
    guinier_grp = _new_process(nxs, parent_grp, name, sequence_index, "freesas.autorg")
    guinier_autorg = nxs.new_class(guinier_grp, "autorg", "NXcollection")
    guinier_gpa = nxs.new_class(guinier_grp, "gpa", "NXcollection")
    guinier_guinier = nxs.new_class(guinier_grp, "guinier", "NXcollection")
    guinier_data = nxs.new_class(guinier_grp, "result", "NXdata")
    guinier_data.attrs["SILX_style"] = NORMAL_STYLE
    guinier_data.attrs["title"] = "Guinier analysis"

    q, I, err = sasm.T[:3]
    try:
        gpa = auto_gpa(sasm)
    except Exception as error:
        guinier_gpa["Failed"] = f"{error.__class__.__name__}: {error}"
        gpa = None
    else:
        _store_rg(guinier_gpa, gpa, radius_unit)

    try:
        guinier = auto_guinier(sasm)
    except Exception as error:
        guinier_guinier["Failed"] = f"{error.__class__.__name__}: {error}"
        guinier = None
    else:
        _store_rg(guinier_guinier, guinier, radius_unit, q)

    try:
        autorg = autoRg(sasm)
    except Exception as error:
        guinier_autorg["Failed"] = f"{error.__class__.__name__}: {error}"
        autorg = None
    else:
        if autorg.Rg < 0:
            guinier_autorg["Failed"] = "No Guinier region found with this algorithm"
            autorg = None
        else:
            _store_rg(guinier_autorg, autorg, radius_unit, q)

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

    # Guinier plot generation:
    mask = (I > 0) & numpy.isfinite(I) & (q > 0) & numpy.isfinite(q)
    if err is not None:
        mask &= (err > 0.0) & numpy.isfinite(err)
    mask = mask.astype(bool)
    if guinier:
        intercept = numpy.log(guinier.I0)
        slope = -(guinier.Rg**2) / 3.0
        invalid = numpy.where(q > 1.5 / guinier.Rg)[0]
        if invalid.size:
            end = invalid[0]
            mask[end:] = False

    q2 = q[mask] ** 2
    logI = numpy.log(I[mask])
    dlogI = abs(err[mask] / I[mask])
    q2_ds = guinier_data.create_dataset("q2", data=q2.astype(numpy.float32))
    q2_ds.attrs["unit"] = radius_unit + "⁻²"
    q2_ds.attrs["long_name"] = f"q² ({radius_unit}⁻²)"
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
    guinier_grp.attrs["default"] = posixpath.relpath(guinier_data.name, guinier_grp.name)
    return guinier


def kratky_plot(nxs, parent_grp, name, sequence_index, sasm, guinier):
    """Build the dimensionless Kratky plot

    :param nxs: opened Nexus file object
    :param parent_grp: group in which the NXprocess is created
    :param name: name of the NXprocess group
    :param sequence_index: index of the process
    :param sasm: 2D array with q, I, sigma as columns
    :param guinier: Guinier fit (RG_RESULT)
    """
    kratky_grp = _new_process(nxs, parent_grp, name, sequence_index, "freesas.autorg")
    kratky_data = nxs.new_class(kratky_grp, "result", "NXdata")
    kratky_data.attrs["SILX_style"] = NORMAL_STYLE
    kratky_data.attrs["title"] = "Dimensionless Kratky plots"
    kratky_grp.attrs["default"] = posixpath.relpath(kratky_data.name, kratky_grp.name)

    q, I, err = sasm.T[:3]
    Rg = guinier.Rg
    I0 = guinier.I0
    xdata = q * Rg
    ydata = xdata * xdata * I / I0
    dy = xdata * xdata * err / I0
    qRg_ds = kratky_data.create_dataset("qRg", data=xdata.astype(numpy.float32))
    qRg_ds.attrs["interpretation"] = "spectrum"
    qRg_ds.attrs["long_name"] = "q·Rg (unit-less)"
    # Nota the "÷" hereafter is the division sign and not the usual slash
    k_ds = kratky_data.create_dataset("q2Rg2I÷I0", data=ydata.astype(numpy.float32))
    k_ds.attrs["interpretation"] = "spectrum"
    k_ds.attrs["long_name"] = "q²Rg²I(q)/I₀"
    ke_ds = kratky_data.create_dataset("errors", data=dy.astype(numpy.float32))
    ke_ds.attrs["interpretation"] = "spectrum"
    kratky_data_attrs = kratky_data.attrs
    kratky_data_attrs["signal"] = posixpath.basename(k_ds.name)
    kratky_data_attrs["axes"] = posixpath.basename(qRg_ds.name)


def invariants(nxs, parent_grp, name, sequence_index, sasm, guinier):
    """Calculate the Rambo-Tainer invariants and the Porod volume

    :param nxs: opened Nexus file object
    :param parent_grp: group in which the NXprocess is created
    :param name: name of the NXprocess group
    :param sequence_index: index of the process
    :param sasm: 2D array with q, I, sigma as columns
    :param guinier: Guinier fit (RG_RESULT)
    :return: 2-tuple with the Rambo-Tainer invariants and the Porod volume
    """
    rti_grp = nxs.new_class(parent_grp, name, "NXprocess")
    rti_grp["sequence_index"] = sequence_index
    rti_grp["program"] = "freesas.invariants"
    rti_grp["version"] = freesas.version
    rti_data = nxs.new_class(rti_grp, "result", "NXdata")

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

    # Porod
    volume = freesas.invariants.calc_Porod(sasm, guinier)
    volume_ds = rti_data.create_dataset("volume", data=volume)
    volume_ds.attrs["unit"] = "nm³"
    volume_ds.attrs["formula"] = "Porod: V = 2*π²I₀²/(sum_q I(q)q² dq)"
    return rti, volume


def bift_analysis(nxs, parent_grp, name, sequence_index, sasm, guinier, radius_unit, curve_data=None):
    """Pair distribution function by Bayesian indirect Fourier transformation, the equivalent of datgnom

    :param nxs: opened Nexus file object
    :param parent_grp: group in which the NXprocess is created
    :param name: name of the NXprocess group
    :param sequence_index: index of the process
    :param sasm: 2D array with q, I, sigma as columns
    :param guinier: Guinier fit (RG_RESULT)
    :param radius_unit: "nm" or "Å"
    :param curve_data: NXdata with the scattering curve, the BIFT fit is overlaid on top of it
    :return: BIFT statistics or None if the IFT failed
    """
    bift_grp = _new_process(nxs, parent_grp, name, sequence_index, "freesas.bift")
    bift_data = nxs.new_class(bift_grp, "result", "NXdata")
    bift_data.attrs["SILX_style"] = NORMAL_STYLE
    bift_data.attrs["title"] = "Pair distance distribution function p(r)"

    cfg_grp = nxs.new_class(bift_grp, "configuration", "NXcollection")
    q, I, err = sasm.T[:3]
    try:
        bo = BIFT(q, I, err)
        cfg_grp["Rg"] = guinier.Rg
        # Pretty limited quality as we have real time constrains
        cfg_grp["npt"] = npt = 64
        cfg_grp["Dmax÷Rg"] = 3
        Dmax = bo.set_Guinier(guinier, Dmax_over_Rg=3)

        # First scan on alpha:
        cfg_grp["alpha_sup"] = alpha_max = bo.guess_alpha_max(npt)
        cfg_grp["alpha_inf"] = 1 / alpha_max
        cfg_grp["alpha_scan_steps"] = 11

        key = bo.grid_scan(Dmax, Dmax, 1, 1.0 / alpha_max, alpha_max, 11, npt)
        Dmax, alpha = key[:2]
        # Then scan on Dmax:
        cfg_grp["Dmax_sup"] = guinier.Rg * 4
        cfg_grp["Dmax_inf"] = guinier.Rg * 2
        cfg_grp["Dmax_scan_steps"] = 5
        key = bo.grid_scan(guinier.Rg * 2, guinier.Rg * 4, 5, alpha, alpha, 1, npt)
        Dmax, alpha = key[:2]
        if bo.evidence_cache[key].converged:
            bo.update_wisdom()
            use_wisdom = True
        else:
            use_wisdom = False
        res = minimize(bo.opti_evidence, (Dmax, log(alpha)), args=(npt, use_wisdom), method="powell")
        cfg_grp["Powell_steps"] = res.nfev
        cfg_grp["Monte-Carlo_steps"] = 0
        stats = bo.calc_stats()
    except Exception as error:
        bift_grp["Failed"] = f"{error.__class__.__name__}: {error}"
        return None

    bift_grp["alpha"] = stats.alpha_avg
    bift_grp["alpha_error"] = stats.alpha_std
    bift_grp["Dmax"] = stats.Dmax_avg
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
    r_ds.attrs["long_name"] = f"radius r({radius_unit})"
    p_ds = bift_data.create_dataset("p(r)", data=stats.density_avg.astype(numpy.float32))
    p_ds.attrs["long_name"] = "Pair distance distribution p(r)"
    p_ds.attrs["interpretation"] = "spectrum"
    bift_data["errors"] = stats.density_std
    bift_data.attrs["signal"] = "p(r)"
    bift_data.attrs["axes"] = "r"

    if curve_data is not None:
        r = stats.radius
        T = numpy.outer(q, r / pi)
        T = (4 * pi * (r[-1] - r[0]) / (len(r) - 1)) * numpy.sinc(T)
        bift_ds = curve_data.create_dataset("BIFT", data=T.dot(stats.density_avg).astype(numpy.float32))
        bift_ds.attrs["interpretation"] = "spectrum"
        curve_data.attrs["auxiliary_signals"] = "BIFT"
    bift_grp.attrs["default"] = posixpath.relpath(bift_data.name, bift_grp.name)
    return stats


def saxs_analysis(nxs, parent_grp, sasm, radius_unit, sequence_index, curve_data=None, first_step=None):
    """Perform the complete analysis: Guinier -> Kratky -> invariants -> BIFT.

    The sequence stops after the Guinier analysis if no Guinier region was found.

    :param nxs: opened Nexus file object
    :param parent_grp: group in which the NXprocess are created
    :param sasm: 2D array with q, I, sigma as columns
    :param radius_unit: "nm" or "Å"
    :param sequence_index: callable (i.e. SequenceIndex) providing the index of each process
    :param curve_data: NXdata with the scattering curve, the BIFT fit is overlaid on top of it
    :param first_step: number used in the name of the first group, the following are incremented.
                       When None, the sequence index is used in the name of the group.
    :return: AnalysisResult
    """
    step = None if first_step is None else itertools.count(first_step)

    def next_process(base):
        "return the name and the sequence index of the next process"
        seq = sequence_index()
        return f"{seq if step is None else next(step)}_{base}", seq

    guinier = guinier_analysis(nxs, parent_grp, *next_process("Guinier_analysis"), sasm, radius_unit)
    if guinier is None:
        return AnalysisResult()
    kratky_plot(nxs, parent_grp, *next_process("dimensionless_Kratky_plot"), sasm, guinier)
    rti, volume = invariants(nxs, parent_grp, *next_process("invariants"), sasm, guinier)
    stats = bift_analysis(nxs, parent_grp, *next_process("indirect_Fourier_transformation"),
                          sasm, guinier, radius_unit, curve_data)
    return AnalysisResult(guinier, rti, volume, stats)
