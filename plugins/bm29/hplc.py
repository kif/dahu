"""Data Analysis plugin for BM29: BioSaxs

* HPLC mode: Rebuild the complete chromatogram and perform basic analysis on it.
"""
from __future__ import annotations

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "08/10/2026"
__status__ = "development"
__version__ = "0.6.0"

import copy
import glob
import json

# from dahu.utils import fully_qualified_name
import logging
import math
import os
import posixpath
import re
import time
import zipfile
from typing import NamedTuple

import freesas
import freesas.autorg
import freesas.cormap
import freesas.invariants
import h5py
import matplotlib.pyplot
import numpy
import pyFAI
import pyFAI.integrator.azimuthal
import pyFAI.units
import scipy.ndimage
import scipy.signal
import sklearn
from freesas.app.extract_ascii import write_ascii
from freesas.containers import UVJuice
from freesas.plot import hplc_plot
from pyFAI.containers import ErrorModel, Integrate1dResult
from pyFAI.method_registry import IntegrationMethod
from sklearn.decomposition import NMF
from urllib3.util import parse_url

from dahu.plugin import Plugin

from .analysis import saxs_analysis
from .common import (
    NORMAL_STYLE,
    SAXS_STYLE,
    Ispyb,
    Sample,
    SequenceIndex,
    cmp_float,
    create_nexus_sample,
    str_,
)
from .icat import send_icat
from .integrate import Accumulators
from .ispyb import IspybConnector
from .nexus import Nexus, get_isotime

logger = logging.getLogger("bm29.hplc")
matplotlib.use("Agg")


class NexusJuice(NamedTuple):
    """All information of an integration file"""

    filename: str
    h5path: str
    npt: int
    unit: pyFAI.units.Unit
    idx: numpy.ndarray
    Isum: numpy.ndarray
    q: numpy.ndarray
    I: numpy.ndarray
    sigma: numpy.ndarray
    poni: str
    mask: numpy.ndarray
    energy: float
    polarization: float
    method: tuple
    sample: Sample
    timestamps: numpy.ndarray
    diode: numpy.ndarray | None = None
    ring_current: numpy.ndarray | None = None
    spottiness: numpy.ndarray | None = None
    isotropic: numpy.ndarray | None = None
    accumulators: Accumulators | None = None
    normalization_factor: float | None = None
    count_time: float | None = None
    image_file: str | None = None

    @classmethod
    def read(cls, filename:str):
        """Extract some NexusJuice from a HDF5 file, alternative constructor

        :param filename: name of the file
        :return: NexusJuice instance
        """
        with Nexus(filename, "r") as nxsr:
            entry_name = nxsr.h5.attrs["default"]
            entry_grp = nxsr.h5[entry_name]
            h5path = entry_grp.name
            nxdata_grp = entry_grp[entry_grp.attrs["default"]]
            assert nxdata_grp.name.endswith("hplc")  # we are reading HPLC data
            signal = nxdata_grp.attrs["signal"]
            axis = nxdata_grp.attrs["axes"]
            Isum = nxdata_grp[signal][()]
            idx = nxdata_grp[axis][()]
            integrated = nxdata_grp.parent["result"]
            signal = integrated.attrs["signal"]
            I = integrated[signal][()]
            axes = integrated.attrs["axes"][-1]
            q = integrated[axes][()]
            sigma = integrated["errors"][()]

            npt = len(q)
            unit = pyFAI.units.to_unit(axes + "_" + integrated[axes].attrs["units"])
            integration_grp = nxdata_grp.parent
            poni = str_(integration_grp["configuration/file_name"][()]).strip()
            if not os.path.exists(poni):
                poni = str_(integration_grp["configuration/data"][()]).strip()
            polarization = integration_grp["configuration/polarization_factor"][()]
            method = IntegrationMethod.select_method(
                **json.loads(integration_grp["configuration/integration_method"][()])
            )[0]
            # Unreduced sums of the azimuthal integration, one line per frame.
            # They allow frames to be merged without re-integrating anything and
            # carry both error models. Unavailable in former files.
            if "accumulators" in integration_grp:
                acc_grp = integration_grp["accumulators"]
                accumulators = Accumulators(
                    sum_signal=acc_grp["sum_signal"][()],
                    sum_normalization=acc_grp["sum_normalization"][()],
                    # sum_normalization2 and count are missing in former files,
                    # where frames can only be merged by summing the accumulators
                    sum_normalization2=acc_grp["sum_normalization2"][()]
                    if "sum_normalization2" in acc_grp
                    else None,
                    sum_variance_azimuthal=acc_grp["sum_variance_azimuthal"][()],
                    count=acc_grp["count"][()] if "count" in acc_grp else None,
                    sum_variance_poisson=acc_grp["sum_variance_poisson"][()]
                    if "sum_variance_poisson" in acc_grp
                    else None,
                )
            else:
                accumulators = None
            # Azimuthal heterogeneity, used to spot menisci and parasitic scattering
            if "anisotropy" in integration_grp:
                aniso_grp = integration_grp["anisotropy"]
                spottiness = aniso_grp["spottiness"][()]
                isotropic = aniso_grp["isotropic"][()]
            else:
                spottiness = isotropic = []
            instrument_grp = nxsr.get_class(entry_grp, class_type="NXinstrument")[0]
            # The beam-stop diode is registered as an NXdetector as well:
            # tell them apart by their content rather than by their position.
            detectors = nxsr.get_class(instrument_grp, class_type="NXdetector")
            detector_grp = next(grp for grp in detectors if "pixel_mask" in grp)
            mask = detector_grp["pixel_mask"].attrs["filename"]
            # The 2D frames are an external link to the file written by LImA: it is
            # the only trace of where the acquisition itself was recorded
            link = detector_grp.get("frames", getlink=True)
            image_file = (os.path.normpath(os.path.join(
                              os.path.dirname(os.path.abspath(filename)), link.filename))
                          if isinstance(link, h5py.ExternalLink) else None)
            count_time = (
                detector_grp["count_time"][()] if "count_time" in detector_grp else None
            )
            diode_grp = next(
                (grp for grp in detectors if "normalization_factor" in grp), None
            )
            normalization_factor = (
                diode_grp["normalization_factor"][()] if diode_grp is not None else None
            )
            source_grp = nxsr.get_class(instrument_grp, class_type="NXsource")[0]
            # Beware: the hard link in the measurement group is spelled "ring_curent"
            ring_current = (
                source_grp["current"][()] if "current" in source_grp else []
            )
            mono_grp = nxsr.get_class(instrument_grp, class_type="NXmonochromator")[0]
            energy = mono_grp["energy"][()]
            #             img_grp = nxsr.get_class(entry_grp["3_time_average"], class_type="NXdata")[0]
            #             image2d = img_grp["intensity_normed"][()]
            #             error2d = img_grp["intensity_std"][()]
            # Read the sample description:
            sample_grp = nxsr.get_class(entry_grp, class_type="NXsample")[0]
            sample_name = posixpath.split(sample_grp.name)[-1]

            buffer = str_(sample_grp["buffer"][()]) if "buffer" in sample_grp else ""
            concentration = (
                sample_grp["concentration"][()] if "concentration" in sample_grp else ""
            )
            description = (
                str_(sample_grp["description"][()]) if "description" in sample_grp else ""
            )
            hplc = sample_grp["hplc"][()] if "hplc" in sample_grp else ""
            temperature = (
                sample_grp["temperature"][()] if "temperature" in sample_grp else ""
            )
            temperature_env = (
                sample_grp["temperature_env"][()]
                if "temperature_env" in sample_grp
                else ""
            )
            sample = Sample(
                sample_name,
                description,
                buffer,
                concentration,
                hplc,
                temperature_env,
                temperature,
            )
            meas_grp = nxsr.get_class(entry_grp, class_type="NXdata")[0]
            timestamps = []
            for ts_name in ("timestamps", "time-stamps"):
                if ts_name in meas_grp:
                    timestamps = meas_grp[ts_name][()]
                    break
            if "diode" in meas_grp:
                diode = meas_grp["diode"][()]
            else:
                diode = []

        return cls( filename=filename,
                    h5path=h5path,
                    npt=npt,
                    unit=unit,
                    idx=idx,
                    Isum=Isum,
                    q=q,
                    I=I,
                    sigma=sigma,
                    poni=poni,
                    mask=mask,
                    energy=energy,
                    polarization=polarization,
                    method=method,
                    sample=sample,
                    timestamps=timestamps,
                    diode=diode,
                    ring_current=ring_current,
                    spottiness=spottiness,
                    isotropic=isotropic,
                    accumulators=accumulators,
                    normalization_factor=normalization_factor,
                    count_time=count_time,
                    image_file=image_file,
                    )

    @classmethod
    def concatenate(cls, juices):
        """Merge the juices of several integration files into a single one, with every
        per-frame series put back in acquisition order, alternative constructor

        Each file only covers a slice of the acquisition and nothing guarantees they are
        read in order, so the frame index is the only trustworthy ordering. Series are
        scattered into arrays of `max(frame_id) + 1` entries: a frame no file provides is
        left at zero, and `idx` lists those actually filled.

        :param juices: iterable of NexusJuice sharing the same radial axis
        :return: a single NexusJuice
        """
        juices = [juice for juice in juices if juice is not None]
        if not juices:
            raise ValueError("No integration file to concatenate")
        first = juices[0]
        for juice in juices[1:]:
            if juice.npt != first.npt or not numpy.allclose(juice.q, first.q):
                raise ValueError(f"{juice.filename} and {first.filename} were integrated "
                                 "on different radial axes, they cannot be concatenated")
        idx = numpy.concatenate([juice.idx for juice in juices])
        order = numpy.argsort(idx, kind="stable")
        idx = numpy.ascontiguousarray(idx[order], dtype=numpy.int64)
        nframes = int(idx[-1]) + 1 if idx.size else 0

        def scatter(chunks, dtype=None):
            "Put the per-frame values of every file back at their frame index"
            if any(chunk is None or len(chunk) == 0 for chunk in chunks):
                return None
            values = numpy.concatenate(chunks)[order]
            out = numpy.zeros((nframes,) + values.shape[1:], dtype=dtype or values.dtype)
            out[idx] = values
            return out

        def series(name, dtype=None):
            "Absent series are handed back empty, as `read` does"
            values = scatter([getattr(juice, name) for juice in juices], dtype)
            return [] if values is None else values

        accumulators = None
        if all(juice.accumulators is not None for juice in juices):
            accumulators = Accumulators(
                **{name: scatter([getattr(juice.accumulators, name) for juice in juices])
                   for name in Accumulators._fields})

        return cls(filename=", ".join(juice.filename for juice in juices),
                   h5path=first.h5path,
                   npt=first.npt,
                   unit=first.unit,
                   idx=idx,
                   Isum=series("Isum", numpy.float64),
                   q=first.q,
                   I=series("I", numpy.float32),
                   sigma=series("sigma", numpy.float32),
                   poni=first.poni,
                   mask=first.mask,
                   energy=first.energy,
                   polarization=first.polarization,
                   method=first.method,
                   sample=first.sample,
                   timestamps=series("timestamps", numpy.float64),
                   diode=series("diode", numpy.float64),
                   ring_current=series("ring_current"),
                   spottiness=series("spottiness"),
                   isotropic=series("isotropic"),
                   accumulators=accumulators,
                   normalization_factor=first.normalization_factor,
                   count_time=first.count_time,
                   image_file=first.image_file,
                   )

    def to_result(self, index, error_model="poisson"):
        """Rebuild a pyFAI Integrate1dResult for one frame, i.e. to merge
        several of them with `union`

        Only `sem` is meaningful on the merged result, never `std`: as soon as the
        curves have been renormalized, `sum_normalization2` is propagated for the
        sake of pyFAI's machinery and no longer carries its statistical meaning.
        For the same reason the azimuthal model, whose crossed term in `union`
        weighs the frames with `sum_normalization2`, is not the default.

        :param index: index of the frame within this file
        :param error_model: "poisson" or "azimuthal", picks which of the two
                            variances accumulated by the integration to use
        :return: Integrate1dResult instance
        """
        acc = self.accumulators
        if acc is None:
            raise ValueError(f"No accumulator stored in {self.filename}")
        model = ErrorModel.parse(error_model)
        sum_variance = (acc.sum_variance_azimuthal if model == ErrorModel.AZIMUTHAL
                        else acc.sum_variance_poisson)
        missing = [name for name, value in (("sum_signal", acc.sum_signal),
                                            ("sum_normalization", acc.sum_normalization),
                                            ("sum_normalization2", acc.sum_normalization2),
                                            (f"sum_variance_{model.name.lower()}", sum_variance),
                                            ("count", acc.count))
                   if value is None]
        if missing:
            raise ValueError(f"Unable to rebuild an Integrate1dResult from {self.filename}: "
                             f"missing {', '.join(missing)}")
        result = Integrate1dResult(numpy.asarray(self.q, dtype=numpy.float64),
                                   numpy.zeros(self.npt, dtype=numpy.float64),
                                   numpy.zeros(self.npt, dtype=numpy.float64))
        result._set_sum_signal(numpy.asarray(acc.sum_signal[index], dtype=numpy.float64))
        result._set_sum_variance(numpy.asarray(sum_variance[index], dtype=numpy.float64))
        result._set_sum_normalization(numpy.asarray(acc.sum_normalization[index], dtype=numpy.float64))
        result._set_sum_normalization2(numpy.asarray(acc.sum_normalization2[index], dtype=numpy.float64))
        result._set_count(numpy.asarray(acc.count[index], dtype=numpy.float64))
        result._set_sem(numpy.zeros(self.npt, dtype=numpy.float64))
        result._set_std(numpy.zeros(self.npt, dtype=numpy.float64))
        result._set_unit(self.unit)
        result._set_polarization_factor(self.polarization)
        result._set_method(self.method)
        result._set_error_model(model)
        return result.__recalculate_means__()


SMOOTHING_ALGORITHMS = ("median", "savgol", "mean", "none")

CHROMATOGRAM_QRANGE = (0.1, 1.0)
"Range of q, in nm⁻¹, summed to build the SAXS chromatogram: best signal to noise"

BACKGROUND_QRANGE = (1.5, 4.0)
"Range of q, in nm⁻¹, where the solvent scatters alone: used to tell solvents apart"

SOLUTE_QRANGE = (0.1, 0.5)
"Range of q, in nm⁻¹, where the solute shows up: used to tell whether anything elutes"

DIODE_FILTER_SIZE = 11
"""Width, in frames, of the filter smoothing the beam-stop diode.

Measured over 177 HPLC runs: the noise of the diode goes into the curves one for one,
and smoothing it away takes some 75 % of the scatter of the chromatogram out. Wider is
not better, as the diode also drifts and steps: beyond ~15 frames the filter starts
cutting into those, and at 31 it biases one run out of six by more than 3σ."""

DIODE_NOISE_LIMIT = 1.0
"Relative noise of the diode, in %, above which normalizing on it makes little sense"


def normalize_chromatogram(signal):
    """Scale a chromatogram between 0 and 1

    Signals of different natures, scattered photons and absorbance at several
    wavelengths, only become comparable on a common plot once rescaled.

    :param signal: the chromatogram as a 1d array
    :return: the same curve, between 0 and 1
    """
    signal = numpy.ascontiguousarray(signal, dtype=numpy.float32)
    low = signal.min()
    span = signal.max() - low
    return (signal - low) / span if span > 0 else numpy.zeros_like(signal)


UV_COUNTER = re.compile(r"^w\d+_(\d+)_\d+$")
"BLISS counter of the UV-Vis spectrometer, the captured group being the wavelength in nm"


def find_bliss_master(image_file):
    """Locate the BLISS master file of an acquisition, starting from one of its images

    LImA writes in a `scan0001` directory sitting next to the master file, which is the
    only HDF5 of its parent directory. The master is not written yet when
    IntegrateMultiframe runs, but it is by the time the chromatogram is built.

    :param image_file: name of one of the HDF5 files written by LImA
    :return: name of the master file, or None when it cannot be pinned down
    """
    if not image_file:
        return None
    scan_dir = os.path.dirname(os.path.abspath(image_file))
    if not os.path.basename(scan_dir).startswith("scan"):
        return None
    dataset_dir = os.path.dirname(scan_dir)
    candidates = glob.glob(os.path.join(dataset_dir, "*.h5"))
    if len(candidates) == 1:
        return candidates[0]
    # Several HDF5 side by side: the master is named after the directory holding it
    expected = os.path.join(dataset_dir, os.path.basename(dataset_dir) + ".h5")
    return expected if expected in candidates else None


def _select_scan(h5file, image_file):
    """Pick the scan of `h5file` which acquired `image_file`

    The detector of a scan is a virtual dataset gathering the LImA files, which is the
    one reliable link between the two.
    """
    scans = [h5file[key] for key in h5file
             if isinstance(h5file[key], h5py.Group) and "measurement" in h5file[key]]
    if image_file:
        target = os.path.basename(image_file)
        for scan in scans:
            for name in scan["measurement"]:
                dataset = scan["measurement"].get(name)
                if isinstance(dataset, h5py.Dataset) and dataset.is_virtual:
                    if any(target == os.path.basename(source.file_name)
                           for source in dataset.virtual_sources()):
                        return scan
    return scans[0] if len(scans) == 1 else None


def read_uv_from_bliss(image_file):
    """Read the UV-Vis chromatograms from the BLISS master file of an acquisition

    Unlike the `.dat` written by the spectrometer, which has a clock of its own, those
    counters are sampled frame per frame on the clock of the scan and need no
    re-alignment.

    The master may still be open for writing by the acquisition, so it is opened
    without locking and everything is guarded: failing to read it is not a reason to
    bring the analysis down.

    :param image_file: name of one of the HDF5 files written by LImA
    :return: UVJuice instance, or None when there is nothing to read
    """
    master_file = find_bliss_master(image_file)
    if master_file is None:
        return None
    try:
        with h5py.File(master_file, "r", locking=False) as h5file:
            scan = _select_scan(h5file, image_file)
            if scan is None:
                logger.info(f"No scan matching {image_file} in {master_file}")
                return None
            measurement = scan["measurement"]
            wavelengths, absorbance = [], []
            for name in sorted(measurement):
                matched = UV_COUNTER.match(name)
                if matched:
                    wavelengths.append(float(matched.group(1)))
                    absorbance.append(measurement[name][()])
            if not wavelengths:
                logger.info(f"No UV-Vis counter in {master_file}{scan.name}")
                return None
            if "elapsed_time" in measurement:
                timestamps = measurement["elapsed_time"][()]
            else:
                epoch = measurement["epoch"][()]
                timestamps = epoch - epoch[0]
            # A scan interrupted half-way leaves counters of uneven length
            size = min(len(timestamps), min(len(i) for i in absorbance))
            return UVJuice(numpy.array(wavelengths),
                           numpy.ascontiguousarray(timestamps[:size]),
                           numpy.vstack([i[:size] for i in absorbance]))
    except Exception as err:
        logger.warning(f"Unable to read the UV-Vis data from {master_file}; "
                       f"{err.__class__.__name__}: {err}")
        return None


def smooth_chromatogram(signal, window, algorithm="median", order=2):
    """Smooth-out a per-frame series: the chromatogram, the diode, ...

    The window is centred and the edges are extended with the closest value, so that
    the first and the last frames are not dragged towards zero.

    :param signal: the chromatogram as 1d array
    :param window: half-width of the window: the filter spans 2*window+1 frames
    :param algorithm: one of SMOOTHING_ALGORITHMS. "median" is insensitive to spikes
                      but leaves steps behind, "savgol" follows the curvature of a
                      drift without eating into it, "mean" is the crudest and "none"
                      hands the signal back untouched
    :param order: order of the polynomial, only used by "savgol"
    :return: the smoothed signal, as a 1d array of float64
    """
    signal = numpy.ascontiguousarray(signal, dtype=numpy.float64)
    width = 2 * int(window) + 1
    if width > signal.size:
        # largest odd window which still fits in the signal
        width = signal.size - 1 + signal.size % 2
    if algorithm == "none" or width < 3:
        return signal
    if algorithm == "median":
        return scipy.ndimage.median_filter(signal, width, mode="nearest")
    if algorithm == "mean":
        return scipy.ndimage.uniform_filter1d(signal, width, mode="nearest")
    if algorithm == "savgol":
        return scipy.signal.savgol_filter(signal, width, min(order, width - 1), mode="nearest")
    raise ValueError(f"Unknown smoothing algorithm {algorithm!r}, expected one of {SMOOTHING_ALGORITHMS}")


def estimate_noise(signal, robust=True):
    """Estimate the standard deviation of the noise of a per-frame series

    Second order differences cancel any linear trend, hence
    E[((sᵢ₋₁ - 2sᵢ + sᵢ₊₁)/√6)²] = σ², whatever the drift of the signal. Unlike the
    residuals to a smoothed curve, this does not depend on the filter which is applied
    afterwards, and does not mix the noise with what that filter fails to follow.

    :param signal: the series as a 1d array
    :param robust: estimate from the median absolute deviation rather than from the
                   root mean square, so that spikes and sharp features do not inflate it
    :return: standard deviation of the noise
    """
    signal = numpy.ascontiguousarray(signal, dtype=numpy.float64)
    if signal.size < 3:
        return 0.0
    residuals = numpy.convolve(signal, [1, -2, 1], "valid") / math.sqrt(6)
    if robust:
        return 1.4826 * numpy.median(abs(residuals - numpy.median(residuals)))
    return residuals.std()


PEAK_PROMINENCE = 10.0
"How far a peak must stand out of the chromatogram, in units of its noise"

PEAK_EDGE_MARGIN = 10
"""Frames at either end of the run where a lone excess is taken for a start-up artefact.

An elution still rising at the last frame is a different matter and is kept: measured
over 174 runs, every spurious flag at an edge was confined to a handful of frames,
while a truncated elution spanned hundreds."""


def contiguous_blocks(indices, minlen=1):
    """Split frame indices into the contiguous stretches they form

    :param indices: sorted 1d array of frame indices
    :param minlen: stretches shorter than this are dropped
    :return: list of (start, stop) slices, stop excluded
    """
    indices = numpy.asarray(indices)
    if indices.size == 0:
        return []
    cuts = numpy.where(numpy.diff(indices) > 1)[0] + 1
    return [(int(block[0]), int(block[-1]) + 1)
            for block in numpy.split(indices, cuts) if len(block) >= minlen]


def search_peaks(intensity, q, wmin=10, prominence=PEAK_PROMINENCE, nsigma=5.0):
    """Label all peak regions of the chromatogram

    Peaks are looked for on the sum over CHROMATOGRAM_QRANGE rather than over the whole
    range: above 1 nm⁻¹ there is solvent and noise and nothing else, and dropping those
    bins buys a factor 3 on the signal to noise. The prominence is given in units of the
    noise of that chromatogram, measured on the data themselves, so the threshold does
    not have to be retuned when the beam or the exposure change.

    Regions are cut at the lowest point between two summits, which keeps them disjoint.
    Taking the width at the base of each peak instead would nest a shoulder inside its
    neighbour, and nested fractions are not fractions.

    A peak still rising at the last frame has no summit, hence no prominence, and
    `find_peaks` is blind to it. `solute_frames` is not, so whatever it flags that no
    peak covers is added as a region of its own, with two exceptions: a flag confined to
    the first or last PEAK_EDGE_MARGIN frames is a start-up artefact rather than an
    elution, and a flag on frames which do not scatter like the solvent of the run is a
    change of solvent. A bubble raises the low to high q ratio just as a solute does, so
    `solute_frames` alone cannot tell them apart; `stationary_frames` can.

    :param intensity: 2D array of shape (nframes, nbins)
    :param q: scattering vector, same unit as CHROMATOGRAM_QRANGE
    :param wmin: minimum width for a peak, smaller ones are discarded
    :param prominence: how far a peak must stand out, in units of the noise
    :param nsigma: width of the gate of `solute_frames`, in robust standard deviations
    :return: (labels, count), as `scipy.ndimage.label` returns
    """
    nframes = intensity.shape[0]
    band = (q >= CHROMATOGRAM_QRANGE[0]) & (q <= CHROMATOGRAM_QRANGE[1])
    chromatogram = intensity[:, band if band.sum() > 1 else slice(None)].sum(axis=-1)
    smooth = smooth_chromatogram(chromatogram, max(wmin // 2, 1))
    noise = estimate_noise(chromatogram)
    regions = []

    peaks, _ = scipy.signal.find_peaks(smooth, prominence=prominence * noise, width=wmin)
    if len(peaks):
        # lowest point between two summits: the natural border between two fractions
        cuts = [int(peak + numpy.argmin(smooth[peak:peaks[i + 1]]))
                for i, peak in enumerate(peaks[:-1])]
        feet = scipy.signal.peak_widths(smooth, peaks, rel_height=0.95)
        for i, peak in enumerate(peaks):
            start = max(int(feet[2][i]), cuts[i - 1] if i else 0)
            stop = min(int(feet[3][i]) + 1, cuts[i] if i + 1 < len(peaks) else nframes)
            if stop - start >= wmin:
                regions.append((start, stop))

    # Safety net: an elution `find_peaks` cannot see for lack of a summit
    stationary = numpy.zeros(nframes, dtype=bool)
    stationary[stationary_frames(intensity, q, nsigma)] = True
    for start, stop in contiguous_blocks(solute_frames(intensity, q, nsigma), wmin):
        if stop <= PEAK_EDGE_MARGIN or start >= nframes - PEAK_EDGE_MARGIN:
            continue
        if stationary[start:stop].mean() < 0.5:
            continue
        if not any(a < stop and start < b for a, b in regions):
            regions.append((start, stop))

    # Labelling by hand rather than through `scipy.ndimage.label`, which would weld
    # two fractions sharing a border into one
    labels = numpy.zeros(nframes, dtype=numpy.int32)
    for index, (start, stop) in enumerate(sorted(regions), start=1):
        labels[start:stop] = index
    return labels, len(regions)

def stationary_frames(intensity, q, nsigma=5.0):
    """Index of the frames whose solvent is the one of the run

    Beyond BACKGROUND_QRANGE the solute scatters almost nothing and the solvent is left
    alone, so the level there identifies the liquid in the capillary. A frame departing
    from the run saw another solvent, a bubble or a glitch, and has nothing to do in a
    background. Being blind to the solute, this tells nothing about the elution.

    :param intensity: 2D array of shape (nframes, nbins)
    :param q: scattering vector, same unit as BACKGROUND_QRANGE
    :param nsigma: width of the gate, in robust standard deviations
    :return: 1d array of frame indices, all of them when the gate cannot be applied
    """
    frames = numpy.arange(intensity.shape[0])
    band = (q >= BACKGROUND_QRANGE[0]) & (q <= BACKGROUND_QRANGE[1])
    if band.sum() < 2:
        logger.warning(f"No q in {BACKGROUND_QRANGE} nm⁻¹: skipping the stationarity gate")
        return frames
    plateau = intensity[:, band].mean(axis=-1)
    median = numpy.median(plateau)
    mad = 1.4826 * numpy.median(abs(plateau - median))
    if mad <= 0:
        return frames
    return frames[abs(plateau - median) < nsigma * mad]


def solute_frames(intensity, q, nsigma=5.0):
    """Index of the frames where something is eluting

    The solute scatters at low q while the solvent is alone at high q, so the ratio of
    the two rises with the solute and does not care about the overall scale. It needs no
    summit either, which is why `search_peaks` falls back on it.

    It covers the whole elution, shoulders included, and that is what a background needs:
    leaving the shoulders of a peak in is worse than leaving the whole peak, because they
    drag the SVD fundamental along without being obvious enough for cormap to reject them.

    :param intensity: 2D array of shape (nframes, nbins)
    :param q: scattering vector, same unit as SOLUTE_QRANGE and BACKGROUND_QRANGE
    :param nsigma: how far above the baseline, in robust standard deviations
    :return: 1d array of frame indices, empty when the test cannot be applied
    """
    low = (q >= SOLUTE_QRANGE[0]) & (q <= SOLUTE_QRANGE[1])
    high = (q >= BACKGROUND_QRANGE[0]) & (q <= BACKGROUND_QRANGE[1])
    if low.sum() < 2 or high.sum() < 2:
        logger.warning(f"No q in {SOLUTE_QRANGE} or {BACKGROUND_QRANGE} nm⁻¹: "
                       "skipping the solute detection")
        return numpy.empty(0, dtype=numpy.int64)
    ratio = intensity[:, low].mean(axis=-1) / intensity[:, high].mean(axis=-1)
    median = numpy.median(ratio)
    mad = 1.4826 * numpy.median(abs(ratio - median))
    if mad <= 0:
        return numpy.empty(0, dtype=numpy.int64)
    return numpy.where(ratio > median + nsigma * mad)[0]


def solute_free_frames(intensity, q, nsigma=5.0):
    """Index of the frames where no solute can be detected at all

    `solute_frames` asks whether something is *clearly* eluting and scales its threshold
    on the spread of the indicator over the whole run. That spread is inflated by the
    elution itself — by a factor 6 in the median over 153 runs, and by several hundred
    when the peak covers much of the acquisition — so "5 deviations" really means some
    30 times the noise. That is the right question for the peak search, which must not
    fire on noise, and the wrong one for a background, which has to be free of any trace
    of solute rather than merely outside the obvious peaks.

    Here the scale is the point to point noise of the indicator, which a smooth elution
    does not inflate, so `nsigma` means what it says.

    :param intensity: 2D array of shape (nframes, nbins)
    :param q: scattering vector, same unit as SOLUTE_QRANGE and BACKGROUND_QRANGE
    :param nsigma: how far above the baseline a frame may sit and still be called clean
    :return: 1d array of frame indices, all of them when the test cannot be applied
    """
    frames = numpy.arange(intensity.shape[0])
    low = (q >= SOLUTE_QRANGE[0]) & (q <= SOLUTE_QRANGE[1])
    high = (q >= BACKGROUND_QRANGE[0]) & (q <= BACKGROUND_QRANGE[1])
    if low.sum() < 2 or high.sum() < 2:
        logger.warning(f"No q in {SOLUTE_QRANGE} or {BACKGROUND_QRANGE} nm⁻¹: "
                       "keeping every frame as a background candidate")
        return frames
    ratio = intensity[:, low].mean(axis=-1) / intensity[:, high].mean(axis=-1)
    noise = estimate_noise(ratio)
    if noise <= 0:
        return frames
    return frames[ratio < numpy.median(ratio) + nsigma * noise]


def build_background(intensity, std=None, keep=0.3, q=None, nsigma=5.0, accumulators=None):
    """
    Build a background from a SVD and search for the frames looking most like the background.

    0. discard the frames whose solvent is not the one of the run (`stationary_frames`)
       and those carrying any trace of solute (`solute_free_frames`)
    1. build a coarse approximation based on the SVD.
    2. measure the distance (cormap) of every single frame to the fundamental of the SVD
    3. average frames that looks most like the coarse approximation (with deviation)

    The two gates of step 0 answer orthogonal questions, and neither replaces cormap:
    one compares levels at high q, where only the solvent speaks, the other the low to
    high q ratio, which the solute raises. A curve offset by a change of solvent keeps
    the shape of a background and cormap lets it through; conversely the shoulders of a
    peak keep the level of a background. Cormap is left to catch what neither sees,
    crystallites or parasitic scattering.

    :param intensity: 2D array of shape (nframes, nbins)
    :param std: same as intensity but with the standard deviation.
    :param keep: fraction of the stationary frames to average (<1!). The gate having
                 already weeded out the odd ones, there is little left to reject here
    :param q: scattering vector; without it the stationarity gate is skipped
    :param nsigma: width of the stationarity gate, in robust standard deviations
    :param accumulators: unreduced sums of the integration, one line per frame. Given
                         those, frames are merged the way pyFAI does rather than with a
                         plain mean of the ratios.
    :return: (bg_avg, bg_std, indexes), each 1d of size nbins. + the index of the frames to keep
    """
    if q is None:
        stationary = numpy.arange(intensity.shape[0])
    else:
        stationary = numpy.intersect1d(stationary_frames(intensity, q, nsigma),
                                       solute_free_frames(intensity, q, nsigma))
    U, S, V = numpy.linalg.svd(intensity[stationary].T, full_matrices=False)
    bg1 = numpy.median(V[0]) * S[0] * U[:, 0]
    Pscore = [
        freesas.cormap.measure_longest(
            numpy.ascontiguousarray(bg1 - i, dtype=numpy.float64)
        )
        for i in intensity[stationary]
    ]
    orderd = numpy.argsort(Pscore)
    nkeep = math.ceil(keep * len(stationary))
    to_keep = numpy.sort(stationary[orderd[:nkeep]])
    if accumulators is None:
        bg_avg = intensity[to_keep].mean(axis=0)
        if std is not None:
            bg_std = numpy.sqrt(((std[to_keep]) ** 2).sum(axis=0)) / len(to_keep)
        else:
            bg_std = None
    else:
        # Summing the unreduced sums is exactly what `Integrate1dResult.union` does for
        # this error model, without the copies: every frame is weighted by its own
        # normalization, which a mean of the ratios silently takes as equal, and the
        # uncertainty comes out as the sem of that weighted mean, diode noise included.
        variance = (accumulators.sum_variance_azimuthal
                    if accumulators.sum_variance_poisson is None
                    else accumulators.sum_variance_poisson)
        normalization = accumulators.sum_normalization[to_keep].sum(axis=0, dtype=numpy.float64)
        bg_avg = accumulators.sum_signal[to_keep].sum(axis=0, dtype=numpy.float64) / normalization
        bg_std = numpy.sqrt(variance[to_keep].sum(axis=0, dtype=numpy.float64)) / normalization
    return bg_avg, bg_std, to_keep, stationary


ZIP_COMPRESSION = zipfile.ZIP_DEFLATED
"""Compression of the archive of curves.

The archive used to be stored uncompressed, which costs some 24 MB for a run of 500
frames; deflate brings it down to 7. bzip2 and lzma do better still, 5.9 and 4.4 MB,
but neither is read by the archive manager shipped with Windows, and these files are
handed to users."""

ZIP_COMPRESSLEVEL = 9
"Deflate level: the slowest, 2 s for a run of 500 frames, which is lost in the noise here"


def save_zip(filename, config, intensity, sigma, background=None, fractions=None,
             dat_template=None):
    """Save the curves of a run into a zip archive

    The archive holds, at its root, the background and one file per fraction, and in a
    `frames` directory the curve of every single frame as it was integrated. Only the
    fractions are background subtracted, the frames are not.

    :param filename: name of the zip-file
    :param config: this is some NexusJuice namedtuple. we use only q and the sample description.
    :param intensity: 2D array with the intensity of the stack of curves
    :param sigma: 2D array with the uncertainties of the stack of frames
    :param background: (intensity, sigma) of the averaged background, or None
    :param fractions: iterable of (first_frame, last_frame, intensity, sigma), subtracted
    :param dat_template: template for the zipped filenames: by default "{basename(filename)}_%04i.dat"
    :return: nothing
    """
    if dat_template is None:
        basename = os.path.basename(filename)
        base = os.path.splitext(basename)[0]
        dat_template =  base + "_%04i.dat"
    common = {"q": config.q}
    if config.sample:
        sample = config.sample
        if sample.name:
            common["sample"] = sample.name
        if sample.buffer:
            common["buffer"] = sample.buffer
        if sample.temperature_env:
            common["storage temperature"] = sample.temperature_env
        if sample.temperature:
            common["exposure temperature"] = sample.temperature
        if sample.concentration:
            common["concentration"] = sample.concentration

    def curve(I, std):
        "A single curve, with the description of the sample attached"
        return dict(common, I=I, std=std)

    with zipfile.ZipFile(filename, "w", compression=ZIP_COMPRESSION,
                         compresslevel=ZIP_COMPRESSLEVEL) as z:
        if background is not None:
            z.writestr("buffer.dat", write_ascii(curve(*background)))
        for first, last, I, std in fractions or ():
            z.writestr(f"fraction_{first}-{last}.dat", write_ascii(curve(I, std)))
        for idx, (i, s) in enumerate(zip(intensity, sigma)):
            z.writestr(posixpath.join("frames", dat_template % idx),
                       write_ascii(curve(i, s)))


class HPLC(Plugin):
    """Rebuild the complete chromatogram and perform basic analysis on it.

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
      "nmf_components": 5,
      "diode_medfilt": 11,        # width of the filter, in frames; 0 or 1 to disable
      "diode_filter": "median",   # or savgol, mean, none
      "uv_datafile": "path to UV .dat file in some gallery",
      "uv_offset": 0.0,           # seconds to add to the UV time-stamps of the .dat
      "background_keep": 0.3,     # fraction of the clean frames averaged as background
      "wait_for": [jobid_img001, jobid_img002],
      "plugin_name": "bm29.hplc"
    }
    """

    NMF_COMP = 5
    "Default number of Non-negative matrix factorization components. Correspond to the number of spices"

    def __init__(self):
        Plugin.__init__(self)
        self.input_files = []
        self.nxs = None
        self.output_file = None
        self.juices = []
        self.juice = None
        self.accumulators = None
        self.uv_source = ""
        self.uv_data = None
        self.nmf_components = self.NMF_COMP
        self.to_pyarch = {}
        self.ispyb = None
        self.sequence_index = SequenceIndex(0)
        self._time_digits = 0

    def setup(self):
        Plugin.setup(self)

        for job_id in self.input.get("wait_for", []):
            self.wait_for(job_id)

        self.input_files = [
            os.path.abspath(i) for i in self.input.get("integrated_files", "")
        ]

        self.output_file = self.input.get("output_file")
        if not self.output_file:
            dirname, basename = os.path.split(
                os.path.commonprefix(self.input_files) + "_hplc.h5"
            )
            dirname = os.path.dirname(dirname)
            #            dirname = os.path.join(dirname, "processed")
            dirname = os.path.join(dirname, "hplc")
            self.output_file = os.path.join(dirname, basename)
            if not os.path.isdir(dirname):
                try:
                    os.makedirs(dirname)
                except Exception as err:
                    self.log_warning(
                        f"Unable to create dir {dirname}. {type(err)}: {err}"
                    )

            self.log_warning("No output file provided, using " + self.output_file)
        self.nmf_components = int(self.input.get("nmf_components", self.NMF_COMP))

        uv_datafile = self.input.get("uv_datafile")
        if uv_datafile and os.path.exists(uv_datafile):
            try:
                self.uv_data = UVJuice.from_file(uv_datafile)
                self.uv_source = uv_datafile
            except Exception as err:
                self.uv_data = None
                self.log_warning(
                    f"Unable to parse {uv_datafile}; {err.__class__.__name__}: {err}"
                )

        # Manage gallery here
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

    def process(self):
        self.create_nexus()
        self.to_pyarch["hdf5_filename"] = self.output_file
        self.to_pyarch["chunk_size"] = self.juices[0].Isum.size
        self.to_pyarch["id"] = os.path.commonprefix(self.input_files)
        self.to_pyarch["sample_name"] = self.juices[0].sample.name
        self.build_plot()
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
        entry_grp = nxs.new_entry(
            "entry",
            self.input.get("plugin_name", "dahu"),
            title="BioSaxs HPLC experiment",
            force_time=get_isotime(),
        )
        entry_grp["version"] = __version__
        nxs.h5.attrs["default"] = entry_grp.name.strip("/")

        # Configuration
        cfg_grp = nxs.new_class(entry_grp, "configuration", "NXnote")
        cfg_grp.create_dataset(
            "data", data=json.dumps(self.input, indent=2, separators=(",\r\n", ":\t"))
        )
        cfg_grp.create_dataset("format", data="text/json")

        # Process 0: Measurement group
        input_grp = nxs.new_class(entry_grp, "0_measurement", "NXcollection")
        input_grp["sequence_index"] = self.sequence_index()

        for idx, filename in enumerate(self.input_files):
            juice = NexusJuice.read(filename)
            if juice is not None:
                rel_path = os.path.relpath(
                    os.path.abspath(filename),
                    os.path.dirname(os.path.abspath(self.output_file)),
                )
                input_grp[f"LImA_{idx:04d}"] = h5py.ExternalLink(rel_path, juice.h5path)
                self.juices.append(juice)

        # Every file covers a slice of the acquisition: put the series back in order
        self.juice = juice = NexusJuice.concatenate(self.juices)

        # The BLISS master file was not written yet when the frames were integrated,
        # but it is by now, and its UV-Vis counters share the clock of the frames.
        # Prefer them over the `.dat`, whose time base has to be realigned by hand.
        from_bliss = read_uv_from_bliss(juice.image_file)
        if from_bliss is not None:
            if self.uv_data:
                self.log_warning(f"Using the UV-Vis counters of the BLISS master file "
                                 f"rather than {self.uv_source}: same time base as the frames")
            self.uv_data = from_bliss
            self.uv_source = find_bliss_master(juice.image_file)
        q = juice.q
        unit = juice.unit
        radial_unit, unit_name = str_(unit).split("_", 1)

        # Sample: outsourced !
        create_nexus_sample(nxs, entry_grp, juice.sample)

        nframes = int(juice.idx[-1]) + 1
        nbin = q.size

        I = juice.I
        sigma = juice.sigma
        Isum = juice.Isum

        ids = numpy.arange(nframes)
        idx = juice.idx
        timestamps = self.to_pyarch["time"] = juice.timestamps

        if len(timestamps):
            self._time_digits = len(f"{timestamps[-1]:.0f}")
        else:
            self._time_digits = 1

        # Process 1: renormalization of the curves on the smoothed diode.
        # The group is created even when the filter leaves the data alone, so that the
        # sequence index of the following steps does not depend on the options.
        diode_raw = juice.diode
        if len(diode_raw) == 0:
            self.log_error("No beam-stop diode in the integrated files: there is "
                           "nothing to normalize the curves with")
        filter_size = self.input.get("diode_medfilt", DIODE_FILTER_SIZE)
        algorithm = self.input.get("diode_filter", "median")
        if filter_size < 2:
            algorithm = "none"
        nrm_grp = nxs.new_class(entry_grp, "1_renormalize", "NXprocess")
        nrm_grp["sequence_index"] = self.sequence_index()
        nrm_grp["filter_used"] = algorithm
        nrm_grp["filter_size"] = filter_size
        diode_smooth = smooth_chromatogram(diode_raw, filter_size // 2, algorithm)
        noise = estimate_noise(diode_raw)
        residual = (diode_raw - diode_smooth).std()
        relative_noise = 100.0 * noise / diode_raw.mean()
        if relative_noise > DIODE_NOISE_LIMIT:
            self.log_warning(f"The beam-stop diode reads {relative_noise:.1f} % of noise, "
                             f"well over the {DIODE_NOISE_LIMIT} % above which normalizing "
                             "on it is meaningless: check the beam and the diode")
        # Uncertainty of the value the curves are divided by: the scatter around the
        # smoothed curve when there is one, the noise of a single reading otherwise.
        diode_error = noise if algorithm == "none" else residual
        noise_ds = nrm_grp.create_dataset(
            "noise", data=relative_noise)
        noise_ds.attrs["unit"] = r"%"
        noise_ds.attrs["formula"] = "1.4826·MAD((dᵢ₋₁-2dᵢ+dᵢ₊₁)/√6) ÷ mean(d)"
        noise_ds.attrs["comment"] = "Noise of a single diode reading, filter independent"
        residual_ds = nrm_grp.create_dataset(
            "residual", data=100.0 * residual / diode_raw.mean())
        residual_ds.attrs["unit"] = r"%"
        residual_ds.attrs["comment"] = ("Scatter of the raw diode around the smoothed one. "
                                        "It exceeds `noise` by whatever the filter cannot "
                                        "follow, and is the honest uncertainty on `smooth`")

        diode_data = nxs.new_class(nrm_grp, "diode", "NXdata")
        raw_ds = diode_data.create_dataset("raw", data=diode_raw.astype(numpy.float32))
        raw_ds.attrs["interpretation"] = "spectrum"
        raw_ds.attrs["long_name"] = "Beam-stop diode intensity"
        smooth_ds = diode_data.create_dataset("smooth", data=diode_smooth.astype(numpy.float32))
        smooth_ds.attrs["interpretation"] = "spectrum"
        smooth_ds.attrs["formula"] = ("left untouched" if algorithm == "none" else
                                      f"{algorithm} filter, {2 * (filter_size // 2) + 1} frames wide")
        smooth_err_ds = diode_data.create_dataset(
            "smooth_errors", data=numpy.full(nframes, diode_error, dtype=numpy.float32))
        smooth_err_ds.attrs["interpretation"] = "spectrum"
        smooth_err_ds.attrs["formula"] = "Incertainty on the smoothed diode value"
        smooth_err_ds.attrs["comment"] = ("`noise` when nothing is smoothed, `residual` "
                                          "otherwise, both in absolute units")
        frame_ds = diode_data.create_dataset("frame_idx", data=ids)
        frame_ds.attrs["interpretation"] = "spectrum"
        frame_ds.attrs["long_name"] = "Frame number"
        # Time on the abscissa, frame_idx kept next to it so that one can switch over
        nrm_time_ds = diode_data.create_dataset("timestamps", data=timestamps,
                                                dtype=numpy.float64)
        nrm_time_ds.attrs["interpretation"] = "spectrum"
        nrm_time_ds.attrs["units"] = "s"
        nrm_time_ds.attrs["long_name"] = "Time (s)"
        diode_data.attrs["axes"] = "timestamps"
        diode_data.attrs["signal"] = "raw"
        diode_data.attrs["auxiliary_signals"] = ["smooth"]
        diode_data.attrs["title"] = "Renormalization"

        # Frames no file provided are left at zero by the concatenation
        scale = numpy.divide(diode_raw, diode_smooth,
                             out=numpy.ones_like(diode_smooth),
                             where=diode_smooth != 0)
        I *= numpy.atleast_2d(scale).T
        Isum *= scale
        sigma *= numpy.atleast_2d(scale).T
        diode = diode_smooth

        # The diode is the dominant source of frame to frame scatter: the detector
        # pixels are many enough for their own variance to average out. Dividing by
        # `d` turns its uncertainty into var(I) = I²·var_d/d², i.e. an extra
        # sum_signal²·var_d/d² on the unreduced variances, as `integrate.py` does
        # in the sample-changer pathway.
        relative_error = numpy.atleast_2d(
            numpy.divide(diode_error, diode_smooth,
                         out=numpy.zeros_like(diode_smooth),
                         where=diode_smooth != 0)).T
        sigma[...] = numpy.hypot(sigma, I * relative_error)

        nrm_data = nxs.new_class(nrm_grp, "result", "NXdata")
        nrm_data.attrs["title"] = ("Curves renormalized on the smoothed diode"
                                   if algorithm != "none" else "Curves as integrated")
        nrm_int_ds = nrm_data.create_dataset(
            "I", data=numpy.ascontiguousarray(I, dtype=numpy.float32), **cmp_float)
        nrm_int_ds.attrs["interpretation"] = "spectrum"
        nrm_int_ds.attrs["units"] = "arbitrary"
        nrm_int_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"
        nrm_std_ds = nrm_data.create_dataset(
            "errors", data=numpy.ascontiguousarray(sigma, dtype=numpy.float32), **cmp_float)
        nrm_std_ds.attrs["interpretation"] = "spectrum"
        nrm_q_ds = nrm_data.create_dataset("q", data=q)
        nrm_q_ds.attrs["interpretation"] = "spectrum"
        nrm_q_ds.attrs["unit"] = unit_name
        nrm_q_ds.attrs["long_name"] = "Scattering vector q (nm⁻¹)"
        nrm_data["timestamps"] = nrm_time_ds
        nrm_data.attrs["signal"] = "I"
        nrm_data.attrs["axes"] = ["timestamps", "q"]
        nrm_data.attrs["SILX_style"] = SAXS_STYLE
        nrm_grp.attrs["default"] = posixpath.relpath(nrm_data.name, nrm_grp.name)

        # Normalizing on the smoothed diode instead of the raw one amounts to
        # scaling the normalization by the inverse factor, like pyFAI's
        # `Integrate1dResult.renormalize` does. The signal, its variance and the
        # pixel count are untouched.
        accumulators = None
        if juice.accumulators is not None:
            raw = juice.accumulators
            inverse = numpy.atleast_2d(numpy.reciprocal(scale)).T
            extra_variance = (raw.sum_signal * relative_error) ** 2
            self.accumulators = accumulators = Accumulators(
                sum_signal=raw.sum_signal,
                sum_normalization=raw.sum_normalization * inverse,
                sum_normalization2=(None if raw.sum_normalization2 is None
                                    else raw.sum_normalization2 * inverse ** 2),
                sum_variance_azimuthal=raw.sum_variance_azimuthal + extra_variance,
                count=raw.count,
                sum_variance_poisson=(None if raw.sum_variance_poisson is None
                                      else raw.sum_variance_poisson + extra_variance))
            acc_grp = nxs.new_class(nrm_grp, "accumulators", "NXcollection")
            acc_grp.attrs["comment"] = (
                "Unreduced sums of the azimuthal integration, one line per frame, "
                "corrected for the renormalization. The intensity of a set of frames "
                "is obtained without re-integrating anything: "
                "sum_signal.sum(axis=0)/sum_normalization.sum(axis=0), or by rebuilding "
                "Integrate1dResult objects and merging them with `union`. Mind that "
                "sum_normalization2 is only propagated to keep pyFAI's machinery happy: "
                "once renormalized it no longer carries its statistical meaning, so use "
                "`sem` and never `std`.")
            acc_grp["q"] = nrm_q_ds
            acc_grp["frame_idx"] = frame_ds
            long_names = {
                "sum_signal": "Σᵢ signalᵢ",
                "sum_normalization": "Σᵢ normalizationᵢ, rescaled on the smoothed diode",
                "sum_normalization2": "Σᵢ normalizationᵢ², rescaled on the smoothed diode",
                "sum_variance_azimuthal": "Σᵢ varianceᵢ, azimuthal error model + diode noise",
                "sum_variance_poisson": "Σᵢ varianceᵢ, poissonian error model + diode noise",
                "count": "Σᵢ pixel countᵢ",
            }
            for name, long_name in long_names.items():
                data = getattr(accumulators, name)
                if data is None:
                    continue
                acc_ds = acc_grp.create_dataset(
                    name, data=numpy.ascontiguousarray(data, dtype=numpy.float32), **cmp_float)
                acc_ds.attrs["interpretation"] = "spectrum"
                acc_ds.attrs["long_name"] = long_name

        # Process 2: Chromatogram
        chroma_grp = nxs.new_class(entry_grp, "2_chromatogram", "NXprocess")
        chroma_grp["sequence_index"] = self.sequence_index()

        # UV-chromatogram
        if self.uv_data:
            uv_data = nxs.new_class(chroma_grp, "UV-Vis", "NXdata")
            uv_data.attrs["title"] = "UV-Vis - Chromatogram"
            uv_data["sequence_index"] = self.sequence_index()
            absorbance = uv_data.create_dataset(
                "absorbance", data=self.uv_data.absorbance
            )
            absorbance.attrs["unit"] = "∅"
            absorbance.attrs["long_name"] = "Absorbance (mAU)"
            absorbance.attrs["interpretation"] = "spectrum"
            absorbance.attrs["SILX_style"] = NORMAL_STYLE
            uv_data.create_dataset("timestamps", data=self.uv_data.timestamps).attrs[
                "unit"
            ] = "s"
            uv_data.create_dataset("wavelengths", data=self.uv_data.wavelengths).attrs[
                "unit"
            ] = "nm"
            uv_data.attrs["signal"] = "absorbance"
            uv_data.attrs["axes"] = ["wavelengths", "timestamps"]

        # SAXS-chromatogram
        hplc_data = nxs.new_class(chroma_grp, "SAXS", "NXdata")
        hplc_data.attrs["title"] = "SAXS - Chromatogram"
        hplc_data["sequence_index"] = self.sequence_index()

        qmin, qmax = CHROMATOGRAM_QRANGE
        band = (q >= qmin) & (q <= qmax)
        band_ds = hplc_data.create_dataset(
            "sum_q_range", data=I[:, band if band.sum() > 1 else slice(None)].sum(axis=-1),
            dtype=numpy.float32)
        band_ds.attrs["interpretation"] = "spectrum"
        band_ds.attrs["long_name"] = f"Summed intensity in the q-range {qmin}-{qmax} nm⁻¹"
        band_ds.attrs["SILX_style"] = NORMAL_STYLE
        band_ds.attrs["comment"] = ("The bins above the range carry solvent and noise only: "
                                    "leaving them out buys a factor 3 on the signal to noise, "
                                    "and this is what the peak search works on")

        sum_ds = hplc_data.create_dataset("sum", data=Isum, dtype=numpy.float32)
        sum_ds.attrs["interpretation"] = "spectrum"
        sum_ds.attrs["long_name"] = "Summed intensity over the whole q-range"
        sum_ds.attrs["SILX_style"] = NORMAL_STYLE

        diode_ds = hplc_data.create_dataset("diode", data=diode, dtype=numpy.float32)
        diode_ds.attrs["interpretation"] = "spectrum"
        diode_ds.attrs["long_name"] = "Beam-stop diode signal"
        diode_ds.attrs["SILX_style"] = NORMAL_STYLE

        frame_ds = hplc_data.create_dataset("frame_ids", data=ids, dtype=numpy.uint32)
        frame_ds.attrs["interpretation"] = "spectrum"
        frame_ds.attrs["long_name"] = "Frame number"

        # `sum` and `diode` stay in the group but out of `auxiliary_signals`: they are
        # three orders of magnitude apart, and overlaying them on a single axis leaves
        # the diode flat on zero and squashes the chromatogram. Whoever wants them
        # together has them, rescaled, in the `result` group next door.
        hplc_data.attrs["signal"] = "sum_q_range"
        hplc_data.attrs["axes"] = "timestamps"  # "frame_ids"
        time_ds = hplc_data.create_dataset(
            "timestamps", data=timestamps, dtype=numpy.float64
        )
        time_ds.attrs["interpretation"] = "spectrum"
        time_ds.attrs["long_name"] = "Time stamps (s)"

        # SAXS and UV-Vis chromatograms overlaid on the frame time base
        merged_data = nxs.new_class(chroma_grp, "result", "NXdata")
        merged_data.attrs["title"] = "Normalized chromatograms"
        merged_data["sequence_index"] = self.sequence_index()
        qmin, qmax = CHROMATOGRAM_QRANGE
        band = (q >= qmin) & (q <= qmax)
        if not band.any():
            self.log_warning(f"No q in [{qmin}, {qmax}] nm⁻¹, summing the whole range instead")
            band = slice(None)
        saxs_ds = merged_data.create_dataset(
            "SAXS", data=normalize_chromatogram(I[:, band].sum(axis=-1)))
        saxs_ds.attrs["interpretation"] = "spectrum"
        saxs_ds.attrs["long_name"] = f"SAXS, Σ I(q) over q ∈ [{qmin}, {qmax}] nm⁻¹"
        auxiliary_signals = []
        if self.uv_data:
            # The spectrometer has its own clock and its own sampling rate: resample it
            # on the frames, padding with zeros wherever it was not recording.
            uv_offset = self.input.get("uv_offset", 0.0)
            uv_timestamps = numpy.asarray(self.uv_data.timestamps, dtype=numpy.float64) + uv_offset
            for wavelength, absorbance in zip(self.uv_data.wavelengths, self.uv_data.absorbance):
                name = f"UV_{wavelength:.0f}nm"
                # Normalize before resampling, so that the padding stays at 0 instead
                # of landing mid-scale and reading as signal
                resampled = numpy.interp(timestamps, uv_timestamps,
                                         normalize_chromatogram(absorbance),
                                         left=0.0, right=0.0)
                uv_ds = merged_data.create_dataset(
                    name, data=numpy.ascontiguousarray(resampled, dtype=numpy.float32))
                uv_ds.attrs["interpretation"] = "spectrum"
                uv_ds.attrs["long_name"] = f"Absorbance at {wavelength:.0f} nm"
                auxiliary_signals.append(name)
            merged_data["uv_offset"] = uv_offset
            merged_data["uv_offset"].attrs["unit"] = "s"
            merged_data["uv_source"] = str(self.uv_source)
        merged_data["timestamps"] = time_ds
        merged_data.attrs["signal"] = "SAXS"
        if auxiliary_signals:
            merged_data.attrs["auxiliary_signals"] = auxiliary_signals
        merged_data.attrs["axes"] = "timestamps"
        merged_data.attrs["SILX_style"] = NORMAL_STYLE
        chroma_grp.attrs["default"] = posixpath.relpath(merged_data.name, chroma_grp.name)
        entry_grp.attrs["default"] = posixpath.relpath(merged_data.name, entry_grp.name)

        # The I(q) themselves are not repeated here: 1_renormalize/result holds them
        chroma_grp.attrs["title"] = str_(self.juices[0].sample)

        # Process 3: SVD decomposition
        svd_grp = nxs.new_class(entry_grp, "3_SVD", "NXprocess")
        svd_grp["sequence_index"] = self.sequence_index()
        logi = numpy.arcsinh(I.T)
        U, S, V = numpy.linalg.svd(logi, full_matrices=False)

        # Number of Eigenvector to keep:
        svd_grp["Ref"] = "https://arxiv.org/pdf/1305.5870.pdf"
        beta = nframes / nbin if nframes <= nbin else 1.0
        omega = 0.56 * beta**3 - 0.95 * beta**2 + 1.82 * beta + 1.43
        tau = numpy.median(S) * omega
        r = numpy.sum(S > tau)

        # Flip axis with negative signal
        flip = V.max(axis=1) < -V.min(axis=1)
        nflip = numpy.where(flip)
        V[nflip] = -V[nflip]
        U[:, nflip] = -U[:, nflip]

        eigen_data = nxs.new_class(svd_grp, "eigenvectors", "NXdata")
        eigen_ds = eigen_data.create_dataset(
            "U", data=numpy.ascontiguousarray(U.T[:r], dtype=numpy.float32)
        )
        eigen_ds.attrs["interpretation"] = "spectrum"
        eigen_ds.attrs["long_name"] = "Eigenvector of the scattering (arcsinh scale)"
        eigen_data["q"] = nrm_q_ds
        eigen_data.attrs["signal"] = "U"
        eigen_data.attrs["axes"] = [".", "q"]
        eigen_data.attrs["SILX_style"] = SAXS_STYLE

        chroma_data = nxs.new_class(svd_grp, "chromatogram", "NXdata")
        chroma_ds = chroma_data.create_dataset(
            "V", data=numpy.ascontiguousarray(V[:r], dtype=numpy.float32)
        )
        chroma_ds.attrs["interpretation"] = "spectrum"
        chroma_ds.attrs["long_name"] = "Weight of the eigenvector along the elution"
        chroma_data["timestamps"] = time_ds
        chroma_data.attrs["signal"] = "V"
        chroma_data.attrs["axes"] = [".", "timestamps"]
        chroma_data.attrs["SILX_style"] = NORMAL_STYLE

        svd_grp.create_dataset("eigenvalues", data=S[:r], dtype=numpy.float32)
        svd_grp.attrs["default"] = posixpath.relpath(chroma_data.name, svd_grp.name)

        # Process 4: NMF matrix decomposition
        nmf_grp = nxs.new_class(entry_grp, "4_NMF", "NXprocess")
        nmf_grp["sequence_index"] = self.sequence_index()
        nmf_grp["program"] = "sklearn.decomposition.NMF"
        nmf_grp["version"] = sklearn.__version__
        nmf = NMF(n_components=self.nmf_components, init="nndsvd", max_iter=1000)
        try:
            W = nmf.fit_transform(I.T)
        except ValueError as err:
            self.log_warning(f"NMF data decomposition failed with: {err}")
            nmf_grp[err.__class__.__name__] = str_(err)
        else:
            eigen_data = nxs.new_class(nmf_grp, "eigenvectors", "NXdata")
            eigen_ds = eigen_data.create_dataset(
                "W", data=numpy.ascontiguousarray(W.T, dtype=numpy.float32)
            )
            eigen_ds.attrs["interpretation"] = "spectrum"
            eigen_ds.attrs["units"] = "arbitrary"
            eigen_ds.attrs["long_name"] = "Scattering of the component"
            eigen_data["q"] = nrm_q_ds
            eigen_data.attrs["signal"] = "W"
            eigen_data.attrs["axes"] = [".", "q"]
            eigen_data.attrs["SILX_style"] = SAXS_STYLE

            H = nmf.components_
            chroma_data = nxs.new_class(nmf_grp, "chromatogram", "NXdata")
            chroma_ds = chroma_data.create_dataset(
                "H", data=numpy.ascontiguousarray(H, dtype=numpy.float32)
            )
            chroma_ds.attrs["interpretation"] = "spectrum"
            chroma_ds.attrs["long_name"] = "Concentration of the component along the elution"
            chroma_data["timestamps"] = time_ds
            chroma_data.attrs["signal"] = "H"
            chroma_data.attrs["axes"] = [".", "timestamps"]
            chroma_data.attrs["SILX_style"] = NORMAL_STYLE
            nmf_grp.attrs["default"] = posixpath.relpath(chroma_data.name, nmf_grp.name)

        # Process 5: Background estimation
        bg_grp = nxs.new_class(entry_grp, "5_background", "NXprocess")
        bg_grp["sequence_index"] = self.sequence_index()
        bg_grp["keep"] = keep = self.input.get("background_keep", 0.3)
        bg_grp["keep"].attrs["info"] = (
            "Fraction of the stationary curves to be considered as background"
        )
        bg_avg, bg_std, to_keep, stationary = build_background(
            I, sigma, keep=keep, q=q, accumulators=accumulators)
        to_keep = numpy.ascontiguousarray(to_keep, dtype=numpy.int32)
        kept_ds = bg_grp.create_dataset("kept", data=to_keep)
        kept_ds.attrs["info"] = (
            "Index of curves used to calculate the background scattering"
        )
        stationary_ds = bg_grp.create_dataset(
            "stationary", data=numpy.ascontiguousarray(stationary, dtype=numpy.int32))
        stationary_ds.attrs["info"] = (
            f"Index of the curves kept as background candidates: those scattering over "
            f"q ∈ {BACKGROUND_QRANGE} nm⁻¹ like the solvent of the run, and free of any "
            f"solute, i.e. without an excess over q ∈ {SOLUTE_QRANGE} nm⁻¹"
        )
        eluting = numpy.setdiff1d(numpy.arange(len(I)), solute_free_frames(I, q))
        eluting_ds = bg_grp.create_dataset(
            "eluting", data=numpy.ascontiguousarray(eluting, dtype=numpy.int32))
        eluting_ds.attrs["info"] = (
            "Index of the curves where something elutes, shoulders included. They are "
            "excluded from the background, whatever the fractions finally retained"
        )
        self.log_warning(f"Background candidates: {len(stationary)} frames out of {len(I)}; "
                         f"{len(eluting)} carry a solute and "
                         f"{len(I) - len(eluting) - len(stationary)} did not scatter like "
                         "the solvent of the run")
        self.to_pyarch["buffer_frames"] = to_keep
        self.to_pyarch["buffer_I"] = bg_avg
        self.to_pyarch["buffer_Stdev"] = bg_std
        bg_data = nxs.new_class(bg_grp, "result", "NXdata")
        bg_data.attrs["signal"] = "I"
        bg_data.attrs["SILX_style"] = SAXS_STYLE
        bg_data.attrs["axes"] = radial_unit
        bg_ds = bg_data.create_dataset(
            "I", data=numpy.ascontiguousarray(bg_avg, dtype=numpy.float32)
        )
        bg_ds.attrs["interpretation"] = "spectrum"
        bg_ds.attrs["units"] = "arbitrary"
        bg_ds.attrs["long_name"] = "Intensity of the background (absolute, normalized on water)"
        bg_q_ds = bg_data.create_dataset(
            radial_unit, data=numpy.ascontiguousarray(q, dtype=numpy.float32)
        )
        bg_q_ds.attrs["units"] = unit_name
        radius_unit = "nm" if "nm" in unit_name else "Å"
        bg_q_ds.attrs["long_name"] = f"Scattering vector q ({radius_unit}⁻¹)"
        bg_std_ds = bg_data.create_dataset(
            "errors", data=numpy.ascontiguousarray(bg_std, dtype=numpy.float32)
        )
        bg_std_ds.attrs["interpretation"] = "spectrum"
        bg_grp.attrs["default"] = posixpath.relpath(bg_data.name, bg_grp.name)
        I_sub = I - bg_avg
        Istd_sub = numpy.sqrt(sigma**2 + bg_std**2)

        self.to_pyarch["scattering_I"] = I
        self.to_pyarch["scattering_Stdev"] = sigma
        self.to_pyarch["subtracted_I"] = I_sub
        self.to_pyarch["subtracted_Stdev"] = Istd_sub
        self.to_pyarch["sum_I"] = Isum

        # Process 6: fraction of chromatogram analysis
        fraction_grp = nxs.new_class(entry_grp, "6_SEC_fractions", "NXprocess")
        fraction_grp["sequence_index"] = self.sequence_index()
        fraction_grp["minimum_size"] = window = 10

        fractions, nfractions = search_peaks(I, q, window)
        self.to_pyarch["merge_frames"] = numpy.zeros((nfractions, 2), dtype=numpy.int32)
        self.to_pyarch["merge_I"] = numpy.zeros((nfractions, nbin), dtype=numpy.float32)
        self.to_pyarch["merge_Stdev"] = numpy.zeros(
            (nfractions, nbin), dtype=numpy.float32
        )

        if nfractions:
            for i, fraction in enumerate(
                scipy.ndimage.find_objects(fractions, nfractions)
            ):
                self.one_fraction(fraction[0], i, nxs, fraction_grp)

        # Now that the background and the fractions are known, the archive can hold them
        save_zip(os.path.splitext(self.output_file)[0] + ".zip", self.juices[0], I, sigma,
                 background=(bg_avg, bg_std),
                 fractions=[(int(first), int(last), merged, deviation)
                            for (first, last), merged, deviation
                            in zip(self.to_pyarch["merge_frames"],
                                   self.to_pyarch["merge_I"],
                                   self.to_pyarch["merge_Stdev"])])

        # Process 7: All other calculation for ISPyB:
        t = self.build_ispyb_group(nxs, entry_grp)
        self.log_warning(f"Ispyb structure creation took {t:.3f}s")

        # Mind to close the file
        nxs.close()

    def one_fraction(self, fraction, index, nxs, top_grp):
        """
        :param fraction: slice with start and end
        :param index: index of the fraction
        :param nxs: opened Nexus file object
        :param top_grp: top level nexus group to start building into.

        """
        q = self.juices[0].q
        unit = self.juices[0].unit
        sample = self.juices[0].sample
        radial_unit, unit_name = str_(unit).split("_", 1)

        I_sub = self.to_pyarch["subtracted_I"]
        sigma = self.to_pyarch["subtracted_Stdev"]

        time = self.to_pyarch["time"]

        template = f"%0{self._time_digits}.0fs-%0{self._time_digits}.0fs"
        time_slice = template % (time[fraction.start],
                                 time[min(fraction.stop, time.size-1)])
        # time_slice = f"{time[fraction.start]:.0f}s-{time[min(fraction.stop, time.size-1)]:.0f}s"
        f_grp = nxs.new_class(
            top_grp,
            time_slice,
            "NXprocess",
        )
        f_grp["sequence_index"] = self.sequence_index()
        f_grp["first_frame"] = fraction.start
        f_grp["last_frame"] = fraction.stop
        # f_grp["program"] = "dahu.plugins.bm29.hplc"
        # f_grp["version"] = __version__
        f_grp["date"] = get_isotime()

        avg_data = nxs.new_class(f_grp, "1_average", "NXdata")
        avg_data["sequence_index"] = self.sequence_index()
        avg_data.attrs["SILX_style"] = SAXS_STYLE
        avg_data.attrs["title"] = (
            f"{sample.name}, frames {fraction.start}-{fraction.stop} averaged ({time_slice}), buffer subtracted"
        )
        avg_data.attrs["signal"] = "I"
        avg_data.attrs["axes"] = radial_unit
        f_grp.attrs["default"] = posixpath.relpath(avg_data.name, f_grp.name)
        avg_q_ds = avg_data.create_dataset(
            radial_unit, data=numpy.ascontiguousarray(q, dtype=numpy.float32)
        )
        avg_q_ds.attrs["units"] = unit_name
        radius_unit = "nm" if "nm" in unit_name else "Å"
        avg_q_ds.attrs["long_name"] = f"Scattering vector q ({radius_unit}⁻¹)"
        if self.accumulators is None:
            I_frc = I_sub[fraction].mean(axis=0)
            fsig2 = sigma[fraction] ** 2
            sigma_frc = numpy.sqrt(fsig2.sum(axis=0)) / fsig2.shape[0]
        else:
            # Merge the frames of the peak on the unreduced sums, exactly as the
            # background is merged, and take the background out once at the end rather
            # than frame by frame: averaging ratios would weigh frames which did not
            # receive the same flux as if they had, and subtracting first would add the
            # uncertainty of the background as many times as there are frames.
            acc = self.accumulators
            variance = (acc.sum_variance_azimuthal
                        if acc.sum_variance_poisson is None
                        else acc.sum_variance_poisson)
            normalization = acc.sum_normalization[fraction].sum(axis=0, dtype=numpy.float64)
            I_frc = (acc.sum_signal[fraction].sum(axis=0, dtype=numpy.float64)
                     / normalization - self.to_pyarch["buffer_I"])
            sigma_frc = numpy.sqrt(
                variance[fraction].sum(axis=0, dtype=numpy.float64) / normalization ** 2
                + self.to_pyarch["buffer_Stdev"] ** 2)
        ai2_int_ds = avg_data.create_dataset(
            "I", data=numpy.ascontiguousarray(I_frc, dtype=numpy.float32)
        )
        ai2_std_ds = avg_data.create_dataset(
            "errors", data=numpy.ascontiguousarray(sigma_frc, dtype=numpy.float32)
        )

        ai2_int_ds.attrs["interpretation"] = "spectrum"
        ai2_int_ds.attrs["units"] = "arbitrary"
        ai2_int_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"
        #  ai2_int_ds.attrs["uncertainties"] = "errors" #this does not work
        ai2_std_ds.attrs["interpretation"] = "spectrum"
        ai2_int_ds.attrs["units"] = "arbitrary"

        self.to_pyarch["merge_frames"][index, 0] = fraction.start
        self.to_pyarch["merge_frames"][index, 1] = fraction.stop
        self.to_pyarch["merge_I"][index] = I_frc
        self.to_pyarch["merge_Stdev"][index] = sigma_frc

        # Process 2-5: Guinier analysis, Kratky plot, invariants and BIFT
        sasm = numpy.vstack((q, I_frc, sigma_frc)).T
        analysis = saxs_analysis(
            nxs, f_grp, sasm, radius_unit, self.sequence_index,
            curve_data=avg_data, first_step=2,
        )
        if analysis.guinier is None:
            f_grp.attrs["default"] = posixpath.relpath(avg_data.name, f_grp.name)
            self.log_error(
                "No Guinier region found, data of dubious quality", do_raise=False
            )
            return
        self.Vc = analysis.rti.Vc
        self.mass = analysis.rti.mass
        if analysis.bift is not None:
            self.Dmax = analysis.bift.Dmax_avg

    def build_ispyb_group(self, nxs, top_grp):
        """Build the ispyb group inside the HDF5/Nexus file and all associated calculation
        :param nsx: opened nexus file
        :param top_grp: top level group where to start working
        :return: runtime
        """

        def normalize(dataset, dtype=numpy.float32):
            if numpy.dtype(dtype).kind == "f":
                "Deal with Nans"
                mask = numpy.logical_not(numpy.isfinite(dataset))
                dataset = dataset.copy()
                dataset[mask] = 0.0
            return numpy.ascontiguousarray(dataset, dtype=dtype)

        keys = [
            "buffer_frames",
            "buffer_I",
            "buffer_Stdev",
            "Dmax",
            "gnom",
            "I0",
            "I0_Stdev",
            "mass",
            "mass_Stdev",
            "merge_frames",
            "merge_I",
            "merge_Stdev",
            "q",
            "Qr",
            "Qr_Stdev",
            "quality",
            "Rg",
            "Rg_Stdev",
            "scattering_I",
            "scattering_Stdev",
            "subtracted_I",
            "subtracted_Stdev",
            "sum_I",
            "time",
            "total",
            "Vc",
            "Vc_Stdev",
            "volume",
        ]
        keys_extra = {
            "Dmax": "Maximum diameter from IFT, skipped",
            "gnom": "Radius of gyration from IFT, skipped",
            "I0": "Forward scattering from GPA",
            "I0_Stdev": "Uncertainty on forward scattering",
            "mass": "Estimated protein weight from RT analysis",
            "mass_Stdev": "Uncertainty on the mass",
            "Qr": "RT invariant",
            "Qr_Stdev": "Uncertainty on RT invariant",
            "quality": "Quality estimated from GPA",
            "Rg": "radius of gyration obtained from GPA",
            "Rg_Stdev": "uncertainty on Rg",
            "total": "skipped",
            "Vc": "Volume of correlation obtained from RT",
            "Vc_Stdev": "Uncertainty on volume",
            "volume": "Molecular volume obtained from Porrod analysis",
            "sum_I": "Total scattering of the frame",
            "time": "Timestamps",
        }
        start_time = time.perf_counter()
        ispyb_grp = nxs.new_class(top_grp, "7_ISPyB", "NXcollection")
        ispyb_grp["sequence_index"] = self.sequence_index()
        ispyb_grp["start_time"] = get_isotime()

        scattering_I = self.to_pyarch["scattering_I"]
        nframes, _nbin = scattering_I.shape

        q = self.juices[0].q.astype(numpy.float64)
        ds = ispyb_grp.create_dataset("q", data=normalize(q, dtype=numpy.float32))
        ds.attrs["info"] = (
            "Scattering vector length, common for all scattering intensities"
        )

        ds = ispyb_grp.create_dataset(
            "buffer_frames", data=self.to_pyarch["buffer_frames"]
        )
        ds.attrs["info"] = (
            "Index of frames used to calculate the background scattring (buffer frames)"
        )
        ds = ispyb_grp.create_dataset(
            "buffer_I", data=normalize(self.to_pyarch["buffer_I"], dtype=numpy.float32)
        )
        ds.attrs["info"] = "Averaged background scattering signal (buffer)"
        ds = ispyb_grp.create_dataset(
            "buffer_Stdev",
            data=normalize(self.to_pyarch["buffer_Stdev"], dtype=numpy.float32),
        )
        ds.attrs["info"] = "Standard deviation of background scattering signal (buffer)"
        ds = ispyb_grp.create_dataset(
            "merge_frames",
            data=normalize(self.to_pyarch["merge_frames"], dtype=numpy.float32),
        )
        ds.attrs["info"] = "Frames merged for each fraction"
        ds = ispyb_grp.create_dataset(
            "merge_I", data=normalize(self.to_pyarch["merge_I"], dtype=numpy.float32)
        )
        ds.attrs["info"] = "Scattering from merged frames in each fraction"
        ds = ispyb_grp.create_dataset(
            "merge_Stdev",
            data=normalize(self.to_pyarch["merge_Stdev"], dtype=numpy.float32),
        )
        ds.attrs["info"] = (
            "Uncertainties on scattering from merged frames in each fraction"
        )
        ds = ispyb_grp.create_dataset(
            "scattering_I",
            data=normalize(self.to_pyarch["scattering_I"], dtype=numpy.float32),
        )
        ds.attrs["info"] = "Scattering of each individual frame"
        ds = ispyb_grp.create_dataset(
            "scattering_Stdev",
            data=normalize(self.to_pyarch["scattering_Stdev"], dtype=numpy.float32),
        )
        ds.attrs["info"] = "Uncertainties on the scattering of each individual frame"
        ds = ispyb_grp.create_dataset(
            "subtracted_I",
            data=normalize(self.to_pyarch["subtracted_I"], dtype=numpy.float32),
        )
        ds.attrs["info"] = "Background subtracted scattering of each individual frame"
        ds = ispyb_grp.create_dataset(
            "subtracted_Stdev",
            data=normalize(self.to_pyarch["subtracted_Stdev"], dtype=numpy.float32),
        )
        ds.attrs["info"] = (
            "Uncertainties on background subtracted scattering of each individual frame"
        )

        for k1 in keys_extra:
            if k1 not in self.to_pyarch:
                self.to_pyarch[k1] = numpy.zeros(nframes, dtype=numpy.float32)
        I_sub = self.to_pyarch["subtracted_I"].astype(numpy.float64)
        Istd_sub = self.to_pyarch["subtracted_Stdev"].astype(numpy.float64)
        for i in range(nframes):
            sasm = numpy.vstack((q, I_sub[i], Istd_sub[i])).T
            try:
                guinier = freesas.autorg.auto_gpa(sasm)
                try:
                    rti = freesas.invariants.calc_Rambo_Tainer(sasm, guinier)
                except Exception:
                    rti = None
                try:
                    porod = freesas.invariants.calc_Porod(sasm, guinier)
                except Exception:
                    porod = None
            except Exception:
                guinier = rti = porod = None
            if guinier is not None:
                for k, v in zip(
                    ["Rg", "Rg_Stdev", "I0", "I0_Stdev", "quality"],
                    ["Rg", "sigma_Rg", "I0", "sigma_I0", "quality"],
                ):
                    self.to_pyarch[k][i] = guinier.__getattribute__(v)
            if rti is not None:
                for k, v in zip(
                    ["Vc", "Vc_Stdev", "Qr", "Qr_Stdev", "mass", "mass_Stdev"],
                    ["Vc", "sigma_Vc", "Qr", "sigma_Qr", "mass", "sigma_mass"],
                ):
                    self.to_pyarch[k][i] = rti.__getattribute__(v)
            if porod is not None:
                self.to_pyarch["volume"][i] = porod
        for k, v in keys_extra.items():
            ds = ispyb_grp.create_dataset(
                k, data=normalize(self.to_pyarch[k], dtype=numpy.float32)
            )
            ds.attrs["info"] = v

        # create all symbolic links at the top level for Ispyb compatibility
        for k in keys:
            nxs.h5[k] = ispyb_grp[k]
        ispyb_grp["end_time"] = get_isotime()
        return time.perf_counter() - start_time

    def build_plot(self):
        """Create a chromatogram in the gallery"""
        gallery = os.path.join(
            os.path.dirname(os.path.abspath(self.output_file)), "gallery"
        )
        filename = os.path.join(gallery, "chromatogram.png")
        if not os.path.isdir(gallery):
            try:
                os.makedirs(gallery)
            except Exception as err:
                self.log_warning(
                    f"Unable to create directory {gallery}; {err.__class__.__name__}: {err}"
                )
                return
        sample = self.to_pyarch.get("sample_name", "sample")
        chromatogram = self.to_pyarch.get("sum_I")

        if chromatogram is not None:
            fractions = self.to_pyarch.get("merge_frames")
            if fractions is not None:
                fractions.sort()
            hplc_plot(
                chromatogram,
                timestamps=self.to_pyarch.get("time"),
                fractions=fractions,
                title=f"Chromatograms of {sample}",
                filename=filename,
                img_format="png",
                uv_data=self.uv_data,
            )
            self.output["chromatogram_file"] = filename

    def send_to_ispyb(self):
        """Data sent to ISPyB are:
        * hdf5File
        * jsonFile built from HDF5
        * hplcPlot various plots generated
        """
        if self.ispyb and self.ispyb.url and parse_url(self.ispyb.url).host:
            ispyb = IspybConnector(*self.ispyb)
            ispyb.send_hplc(self.to_pyarch)
        else:
            self.log_warning(f"Not sending to ISPyB: no valid URL in {self.ispyb}")

    def send_to_icat(self):
        to_icat = copy.copy(self.to_pyarch)
        to_icat["experiment_type"] = "hplc"
        to_icat["sample"] = self.juices[0].sample
        if "volume" in to_icat:
            to_icat.pop("volume")
        metadata = {"scanType": "hplc"}
        gallery = self.ispyb.gallery or os.path.join(
            os.path.dirname(os.path.abspath(self.output_file)), "gallery"
        )
        self.save_csv(
            os.path.join(gallery, "chromatogram.csv"),
            to_icat.get("sum_I"),
            to_icat.get("Rg"),
        )

        if not (self.ispyb.url and parse_url(self.ispyb.url).host):
            self.log_warning("Not sending to iCat: ISPyB metadata not valid")
            return

        return send_icat(
            sample=self.juices[0].sample,
            raw=os.path.dirname(os.path.abspath(self.input_files[0])),
            path=os.path.dirname(os.path.abspath(self.output_file)),
            data=to_icat,
            dataset="HPLC",
            gallery=gallery,
            metadata=metadata,
        )

    def save_csv(self, filename, sum_I, Rg):
        dirname = os.path.dirname(filename)
        if not os.path.isdir(dirname):
            os.makedirs(dirname, exist_ok=True)
        lines = ["id,ΣI,Rg"]
        for idx, (I, rg) in enumerate(zip(sum_I, Rg)):
            lines.append(f"{idx},{I},{rg}")
        lines.append("")
        with open(filename, "w") as csv:
            csv.write(os.linesep.join(lines))
