#!/usr/bin/env python3
#

"""Non-regression tests for the SAXS analysis of BM29: Guinier -> Kratky -> invariants -> BIFT

The HPLC and SubtractBuffer plugins are run on synthetic data (form factor of a sphere)
and the content of the produced HDF5 files is compared to a reference which was produced
with the code prior to the extraction of the analysis into the `analysis` module.

To regenerate the reference (only when the output is expected to change!):

    python -m dahu.plugins.bm29.test.test_analysis --regenerate

The BIFT of the `02s-12s` fraction fits badly on purpose (Chi2r ~ 6e3), so its Powell
descent is sensitive to the version of freesas: expect that fraction, and it alone, to
drift when freesas is upgraded. The reference was last produced with freesas 2026.10.0.
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "08/10/2026"
__status__ = "development"

import json
import logging
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

import numpy

logger = logging.getLogger(__name__)

try:
    import h5py
    import pyFAI.detectors
    import pyFAI.method_registry
    import pyFAI.units
    from pyFAI.integrator.azimuthal import AzimuthalIntegrator

    from .. import analysis, hplc, subtracte
    from ..common import Ispyb, Sample, SequenceIndex
    from ..nexus import Nexus
except ImportError as err:
    logger.warning(f"Unable to import dependencies of BM29 plugins: {err}")
    analysis = hplc = subtracte = None

REFERENCE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference_analysis.json")

# Datasets/attributes which depend on the time, the versions of the libraries or the computer
VOLATILE = ("date", "version", "start_time", "end_time", "file_time", "file_name", "data", "integration_method")
# Messages of exceptions depend on the version of the libraries: only their presence is checked
FAILED = "Failed"


def sphere(q, R=3.0, I0=100.0):
    """Scattering of a sphere of radius R (Rg = R·√(3/5)), with a small background

    :param q: scattering vector in nm⁻¹
    :param R: radius in nm
    :param I0: forward scattering
    :return: intensity
    """
    x = q * R
    return I0 * (3 * (numpy.sin(x) - x * numpy.cos(x)) / x**3) ** 2 + 1e-3


def summarize(filename):
    """Summarize the content of an HDF5 file in a JSON-serializable dict

    :param filename: name of the HDF5 file
    :return: dict with path as key and the summary of the group/dataset as value
    """
    def conv(value):
        if isinstance(value, bytes):
            return value.decode()
        if isinstance(value, numpy.ndarray):
            if value.dtype.kind in "SOU":
                return [conv(i) for i in value.tolist()]
            if value.size == 0:
                return []
            value = value.astype(numpy.float64)
            # Summary statistics of numerical arrays
            return {"shape": list(value.shape),
                    "first": value.ravel()[0], "last": value.ravel()[-1],
                    "mean": value.mean(), "std": value.std(),
                    "min": value.min(), "max": value.max()}
        if isinstance(value, numpy.generic):
            return value.item()
        return value

    res = {}

    def visit(name, obj):
        basename = name.split("/")[-1]
        if basename in VOLATILE:
            return
        entry = {"attrs": {k: conv(obj.attrs[k]) for k in sorted(obj.attrs) if k not in VOLATILE}}
        if isinstance(obj, h5py.Dataset):
            entry["value"] = "" if basename == FAILED else conv(obj[()])
        res[name] = entry

    with h5py.File(filename, "r") as h5:
        h5.visititems(visit)
    return res


def run_hplc(workdir, hplc_module=None, seed=0):
    """Analyse 2 fractions of a synthetic chromatogram: one valid and one without Guinier region

    :param workdir: directory where to write
    :param hplc_module: module containing the HPLC plugin (allows to test an other version)
    :param seed: seed for the random number generator
    :return: summary of the HDF5 file, plugin attributes
    """
    hplc_module = hplc_module or hplc
    rng = numpy.random.default_rng(seed)
    q = numpy.linspace(0.05, 4.0, 500)
    nframes = 20
    intensity = numpy.array([sphere(q) * (1 + 0.01 * rng.standard_normal(q.size)) for _ in range(nframes)])
    sigma = 0.01 * intensity + 1e-4

    plugin = hplc_module.HPLC()
    plugin.juices = [mock.Mock(q=q, unit=pyFAI.units.to_unit("q_nm^-1"), sample=Sample("lysozyme"))]
    plugin.to_pyarch = {"subtracted_I": intensity,
                        "subtracted_Stdev": sigma,
                        "time": numpy.arange(nframes, dtype=numpy.float64),
                        "merge_frames": numpy.zeros((2, 2), dtype=int),
                        "merge_I": numpy.zeros((2, q.size)),
                        "merge_Stdev": numpy.zeros((2, q.size))}
    plugin._time_digits = 2
    filename = os.path.join(workdir, "hplc.h5")
    nxs = Nexus(filename, mode="w")
    entry = nxs.new_entry("entry", "test", force_time="2026-01-01T00:00:00")
    plugin.one_fraction(slice(2, 12), 0, nxs, entry)
    # This fraction has no signal: the Guinier analysis fails
    plugin.to_pyarch["subtracted_I"] = -numpy.ones_like(intensity)
    plugin.one_fraction(slice(13, 18), 1, nxs, entry)
    nxs.close()
    attrs = {k: getattr(plugin, k, None) for k in ("Vc", "mass", "Dmax")}
    return summarize(filename), attrs


def run_subtract(workdir, subtract_module=None, seed=1, background=1000.0, signal=10.0):
    """Run the SubtractBuffer plugin on synthetic data (2 buffers, 1 sample)

    Images are photon-counting like: Poissonian noise, the variance being the number of counts.
    Each buffer has its own realization of the noise.
    Synthetic images are integrated with pyFAI to provide the accumulators; the reading of input files is mocked.

    :param workdir: directory where to write
    :param subtract_module: module containing the SubtractBuffer plugin (allows to test an other version)
    :param seed: seed for the random number generator
    :param background: number of counts per pixel scattered by the buffer
    :param signal: forward scattering of the sample relative to the buffer (negative to have no Guinier region)
    :return: summary of the HDF5 file, plugin attributes, plugin
    """
    subtract_module = subtract_module or subtracte
    rng = numpy.random.default_rng(seed)
    detector = pyFAI.detectors.Detector(pixel1=172e-6, pixel2=172e-6, max_shape=(200, 200))
    ai = AzimuthalIntegrator(dist=1.0, poni1=0.0, poni2=0.0, detector=detector, wavelength=1e-10)
    unit = pyFAI.units.to_unit("q_nm^-1")
    npt = 300
    method = pyFAI.method_registry.IntegrationMethod.select_one_available(("no", "histogram", "cython"), dim=1)
    qimg = ai.array_from_unit(unit=unit)

    def juice(sample, img):
        res = ai._integrate1d_ng(img, npt, variance=numpy.maximum(img, 0), polarization_factor=0.9,
                                 unit=unit, method=method)
        return subtract_module.NexusJuice(filename="dummy.h5", h5path="/entry_0000", npt=npt, unit=unit,
                                          q=res.radial, I=res.intensity, sigma=res.sigma,
                                          poni="dummy.poni", mask="", energy=12.4, polarization=0.9,
                                          method=method, sum_signal=res.sum_signal, sum_variance=res.sum_variance,
                                          sum_normalization=res.sum_normalization, sample=sample,
                                          I_all=[res.intensity], sigma_all=[res.sigma],
                                          sum_normalization2=res.sum_normalization2, count=res.count,
                                          error_model=res.error_model.name)

    def counts(expected):
        "Photon counting: Poissonian noise"
        return rng.poisson(numpy.maximum(expected, 0)).astype(numpy.float64)

    buffer_expected = numpy.full(detector.shape, background)
    form_factor = sphere(qimg, I0=1.0) - 1e-3  # without the constant background of `sphere`
    sample_expected = buffer_expected + background * signal * form_factor
    juices = {"sample.h5": juice(Sample("lysozyme", buffer="tris", concentration=1.0), counts(sample_expected)),
              "buffer_1.h5": juice(Sample("buffer", buffer="tris"), counts(buffer_expected)),
              "buffer_2.h5": juice(Sample("buffer", buffer="tris"), counts(buffer_expected))}
    for name in juices:
        open(os.path.join(workdir, name), "w").close()

    plugin = subtract_module.SubtractBuffer()
    plugin.input = {"plugin_name": "test", "fidelity": 0.0}
    plugin.sample_file = os.path.join(workdir, "sample.h5")
    plugin.buffer_files = [os.path.join(workdir, "buffer_1.h5"), os.path.join(workdir, "buffer_2.h5")]
    plugin.output_file = os.path.join(workdir, "subtracted.h5")
    plugin.ispyb = Ispyb._fromdict({"gallery": workdir})
    plugin.sample_juice = juices["sample.h5"]
    with mock.patch.object(subtract_module.NexusJuice, "read",
                           side_effect=lambda filename: juices[os.path.basename(filename)]):
        try:
            plugin.create_nexus()
        finally:
            plugin.nxs.close()
    attrs = {k: getattr(plugin, k) for k in ("Rg", "I0", "Dmax", "Vc", "mass")}
    attrs["volume"] = plugin.to_pyarch.get("volume")
    attrs["to_pyarch"] = sorted(plugin.to_pyarch)
    return summarize(plugin.output_file), attrs, plugin


def to_json(obj):
    "Convert numpy scalars to plain python objects"
    return json.loads(json.dumps(obj, default=lambda o: o.item() if isinstance(o, numpy.generic) else str(o)))


def generate_reference(filename=REFERENCE, hplc_module=None, subtract_module=None):
    """Generate the reference file, to be used only when the output is expected to change

    :param filename: name of the JSON file
    :param hplc_module, subtract_module: modules to use, by default the current ones
    """
    workdir = tempfile.mkdtemp(prefix="bm29_analysis_")
    try:
        os.makedirs(os.path.join(workdir, "hplc"))
        os.makedirs(os.path.join(workdir, "subtract"))
        hplc_h5, hplc_attrs = run_hplc(os.path.join(workdir, "hplc"), hplc_module)
        sub_h5, sub_attrs, _ = run_subtract(os.path.join(workdir, "subtract"), subtract_module)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    ref = {"hplc": {"hdf5": hplc_h5, "attrs": hplc_attrs},
           "subtract": {"hdf5": sub_h5, "attrs": sub_attrs}}
    with open(filename, "w") as f:
        json.dump(to_json(ref), f, indent=1, sort_keys=True, ensure_ascii=False)


@unittest.skipIf(analysis is None, "Dependencies of BM29 plugins are missing")
class TestAnalysis(unittest.TestCase):
    RTOL = 1e-3
    ATOL = 1e-6

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        with open(REFERENCE) as f:
            cls.reference = json.load(f)

    def setUp(self):
        self.workdir = tempfile.mkdtemp(prefix="bm29_analysis_")

    def tearDown(self):
        shutil.rmtree(self.workdir, ignore_errors=True)

    def assert_close(self, obtained, expected, path):
        "Recursively compare JSON-like structures, with a tolerance on floats"
        if isinstance(expected, dict):
            self.assertIsInstance(obtained, dict, path)
            self.assertEqual(sorted(obtained), sorted(expected), f"keys of {path}")
            for key in expected:
                self.assert_close(obtained[key], expected[key], f"{path}/{key}")
        elif isinstance(expected, list):
            self.assertIsInstance(obtained, list, path)
            self.assertEqual(len(obtained), len(expected), f"length of {path}")
            for idx, (o, e) in enumerate(zip(obtained, expected)):
                self.assert_close(o, e, f"{path}[{idx}]")
        elif isinstance(expected, float) or (isinstance(expected, int) and isinstance(obtained, float)):
            self.assertTrue(numpy.allclose(obtained, expected, rtol=self.RTOL, atol=self.ATOL, equal_nan=True),
                            f"{path}: obtained {obtained}, expected {expected}")
        else:
            self.assertEqual(obtained, expected, path)

    def compare_hdf5(self, obtained, expected):
        missing = sorted(set(expected) - set(obtained))
        extra = sorted(set(obtained) - set(expected))
        self.assertEqual(missing, [], "HDF5 objects missing")
        self.assertEqual(extra, [], "Unexpected HDF5 objects")
        for path in expected:
            self.assert_close(obtained[path], expected[path], path)

    def test_hplc(self):
        "HPLC.one_fraction produces the same output as the reference"
        h5, attrs = run_hplc(self.workdir)
        ref = self.reference["hplc"]
        self.compare_hdf5(h5, ref["hdf5"])
        self.assert_close(to_json(attrs), ref["attrs"], "attrs")

    def test_subtract(self):
        "SubtractBuffer.create_nexus produces the same output as the reference"
        h5, attrs, _ = run_subtract(self.workdir)
        ref = self.reference["subtract"]
        self.compare_hdf5(h5, ref["hdf5"])
        self.assert_close(to_json(attrs), ref["attrs"], "attrs")
        # Physical sanity, independent from the reference: sphere of radius 3 nm, I0 = 10 x buffer
        self.assertAlmostEqual(attrs["Rg"], 3.0 * numpy.sqrt(3 / 5), delta=0.1, msg="Rg of a sphere")
        self.assertAlmostEqual(attrs["I0"], 10000, delta=100, msg="I0")
        self.assertAlmostEqual(attrs["Dmax"], 6.0, delta=0.5, msg="Dmax of a sphere")

    def test_subtract_bift_failure(self):
        "A BIFT which fails is recorded as failed instead of crashing"
        with mock.patch.object(analysis, "BIFT", side_effect=RuntimeError("Simulated failure of the BIFT")):
            h5, attrs, plugin = run_subtract(self.workdir)
        bift = [k for k in h5 if k.endswith("indirect_Fourier_transformation")]
        self.assertEqual(len(bift), 1, "one BIFT group")
        self.assertIn(f"{bift[0]}/{FAILED}", h5, "BIFT failure is recorded")
        self.assertIsNone(attrs["Dmax"], "no Dmax")
        self.assertNotIn("bift", plugin.to_pyarch)
        self.assertAlmostEqual(attrs["Rg"], 3.0 * numpy.sqrt(3 / 5), delta=0.05, msg="Guinier still valid")

    def test_subtract_guinier_failure(self):
        "Without Guinier region, SubtractBuffer raises and the analysis stops after the Guinier step"
        with self.assertRaises(RuntimeError):
            run_subtract(self.workdir, signal=-20.0)
        with h5py.File(os.path.join(self.workdir, "subtracted.h5"), "r") as h5:
            entry = h5[h5.attrs["default"]]
            self.assertTrue(any(k.endswith("Guinier_analysis") for k in entry))
            self.assertFalse(any(k.endswith("Kratky_plot") for k in entry))
            self.assertTrue(entry.attrs["default"].endswith("buffer_subtraction/result"))

    def test_saxs_analysis_naming(self):
        "Groups are named after the sequence index or after a local step number"
        rng = numpy.random.default_rng(2)
        q = numpy.linspace(0.05, 4.0, 500)
        intensity = sphere(q) * (1 + 0.01 * rng.standard_normal(q.size))
        sasm = numpy.vstack((q, intensity, 0.01 * intensity + 1e-4)).T
        with Nexus(os.path.join(self.workdir, "analysis.h5"), mode="w") as nxs:
            entry = nxs.new_entry("entry", "test", force_time="2026-01-01T00:00:00")
            grp1 = nxs.new_class(entry, "global", "NXcollection")
            grp2 = nxs.new_class(entry, "local", "NXcollection")
            res1 = analysis.saxs_analysis(nxs, grp1, sasm, "nm", SequenceIndex(7))
            res2 = analysis.saxs_analysis(nxs, grp2, sasm, "nm", SequenceIndex(7), first_step=2)
            steps = ["Guinier_analysis", "dimensionless_Kratky_plot", "invariants", "indirect_Fourier_transformation"]
            self.assertEqual(set(grp1), {f"{i + 7}_{s}" for i, s in enumerate(steps)})
            self.assertEqual(set(grp2), {f"{i + 2}_{s}" for i, s in enumerate(steps)})
            self.assertEqual([grp2[f"{i + 2}_{s}/sequence_index"][()] for i, s in enumerate(steps)], [7, 8, 9, 10])
        for res in (res1, res2):
            self.assertAlmostEqual(res.guinier.Rg, 3.0 * numpy.sqrt(3 / 5), delta=0.05, msg="Rg of a sphere")
            self.assertAlmostEqual(res.guinier.I0, 100, delta=1, msg="I0")
            self.assertIsNotNone(res.rti)
            self.assertIsNotNone(res.bift)
            # Nota: the fast BIFT (real-time constrains) over-estimates Dmax (2R=6nm) on this curve
            self.assertGreater(res.bift.Dmax_avg, 2 * res.guinier.Rg, msg="Dmax larger than 2·Rg")
        self.assertEqual(res1.guinier.Rg, res2.guinier.Rg)


class TestBackgroundSelection(unittest.TestCase):
    """The scoring and the thresholding of the background candidates in `hplc`.

    Neither is covered by the non-regression reference, which only replays
    `HPLC.one_fraction`, so they are checked here on their own.
    """

    def test_knee_keeps_everything_without_a_shape(self):
        "Degenerate inputs must not throw frames away"
        for score, why in ((numpy.arange(5.0), "fewer values than the floor"),
                           (numpy.ones(50), "no spread at all"),
                           (numpy.r_[numpy.arange(49.0), numpy.nan], "a non finite value")):
            self.assertTrue(hplc.knee_threshold(score).all(), msg=why)

    def test_knee_respects_its_floor(self):
        "A knee falling below BACKGROUND_MIN_FRAMES is overridden by the floor"
        # three good curves and a cliff: the knee alone would keep the three
        score = numpy.r_[0.0, 1.0, 2.0, 1000.0 + numpy.arange(47.0)]
        self.assertEqual(hplc.knee_threshold(score, minimum=1).sum(), 3, "where the knee falls")
        kept = hplc.knee_threshold(score, minimum=8)
        self.assertEqual(kept.sum(), 8, "the floor overrides the knee")
        self.assertTrue(kept[:8].all(), "the floor keeps the lowest scores")

    def test_knee_cuts_at_a_jump(self):
        "A flat population followed by a clearly worse one is split between the two"
        score = numpy.r_[numpy.linspace(1.0, 1.1, 40), numpy.linspace(8.0, 9.0, 10)]
        kept = hplc.knee_threshold(score)
        self.assertEqual(kept.sum(), 40, "the ten outliers are dropped")
        self.assertTrue(kept[:40].all(), "and only those")

    def test_knee_keeps_ties_together(self):
        "Tied scores share a fate: this is what cormap could not do"
        score = numpy.r_[numpy.zeros(30), numpy.full(10, 5.0), numpy.full(10, 50.0)]
        kept = hplc.knee_threshold(score)
        for value in numpy.unique(score):
            self.assertIn(len(set(kept[score == value].tolist())), (1,),
                          msg=f"frames scoring {value} are not treated alike")

    def test_knee_is_order_independent(self):
        "The mask follows the frames, whatever order they come in"
        rng = numpy.random.default_rng(3)
        score = numpy.r_[rng.normal(1.0, 0.05, 60), rng.normal(9.0, 0.5, 15)]
        shuffle = rng.permutation(score.size)
        self.assertTrue(numpy.array_equal(hplc.knee_threshold(score)[shuffle],
                                          hplc.knee_threshold(score[shuffle])))

    def test_zscore_is_null_on_the_reference_itself(self):
        "A curve equal to the reference deviates from it by nothing"
        q = numpy.linspace(0.05, 4.0, 300)
        reference = sphere(q)
        intensity = numpy.vstack([reference, reference])
        std = 0.01 * intensity
        numpy.testing.assert_allclose(hplc.background_zscore(intensity, reference, std), 0.0,
                                      atol=1e-12)

    def test_zscore_counts_sigmas(self):
        "An offset of k sigma scores k, whatever the shape of the curve"
        q = numpy.linspace(0.05, 4.0, 300)
        reference = sphere(q)
        std = 0.01 * reference
        intensity = numpy.vstack([reference + k * std for k in (0.5, 2.0, -3.0)])
        numpy.testing.assert_allclose(
            hplc.background_zscore(intensity, reference, numpy.vstack([std] * 3)),
            [0.5, 2.0, 3.0], rtol=1e-9)

    def test_zscore_sees_a_deviation_that_cancels(self):
        "A tilt of zero net offset must score as high as a uniform one of the same size"
        q = numpy.linspace(0.05, 4.0, 300)
        reference = sphere(q)
        std = 0.01 * reference
        tilt = numpy.ones(q.size); tilt[:q.size // 2] = -1.0     # sums to zero
        intensity = numpy.vstack([reference + 2.0 * std * tilt, reference + 2.0 * std])
        score = hplc.background_zscore(intensity, reference, numpy.vstack([std] * 2))
        self.assertAlmostEqual(score[0], score[1], places=9,
                               msg="the score cancels signs instead of accumulating them")
        self.assertAlmostEqual(score[0], 2.0, places=9)

    def test_zscore_floor_on_pure_noise(self):
        "Scattering alone puts the score at sqrt(2/pi), not at zero"
        rng = numpy.random.default_rng(8)
        reference = numpy.ones(4000)
        std = numpy.full((1, 4000), 0.01)
        intensity = reference + std * rng.standard_normal((1, 4000))
        self.assertAlmostEqual(float(hplc.background_zscore(intensity, reference, std)[0]),
                               numpy.sqrt(2 / numpy.pi), delta=0.03)

    def test_zscore_leans_on_the_reliable_bins(self):
        "The same absolute excess weighs less where the uncertainty is larger"
        reference = numpy.ones(100)
        std = numpy.ones(100)
        std[:50] = 10.0                      # half the range is ten times noisier
        excess = numpy.full((1, 100), 0.1)
        noisy = hplc.background_zscore(reference + excess, reference, std[None, :])
        flat = hplc.background_zscore(reference + excess, reference, numpy.ones((1, 100)))
        self.assertLess(noisy[0], flat[0], "the noisy half must not weigh as much")

    def test_background_rejects_an_obvious_intruder(self):
        "An curve carrying a solute is kept out of the background it would spoil"
        rng = numpy.random.default_rng(4)
        q = numpy.linspace(0.05, 4.0, 300)
        solvent = 100.0 + numpy.zeros_like(q)
        intensity = solvent * (1 + 0.002 * rng.standard_normal((40, q.size)))
        intensity[17] += 20.0 * sphere(q) / sphere(q).max()      # one frame with a solute
        std = 0.002 * intensity
        bg_avg, bg_std, to_keep, candidates, _ = hplc.build_background(intensity, std)
        self.assertNotIn(17, to_keep, "the contaminated frame made it into the background")
        self.assertGreater(len(to_keep), len(intensity) // 2,
                           "a clean pool must not be decimated")
        self.assertEqual(bg_avg.size, q.size)
        self.assertEqual(bg_std.size, q.size)

    def test_background_keep_is_a_ceiling(self):
        "`keep` caps the pool, and only ever downwards"
        rng = numpy.random.default_rng(5)
        intensity = 100.0 * (1 + 0.002 * rng.standard_normal((60, 200)))
        std = 0.002 * intensity
        wide = hplc.build_background(intensity, std)[2]
        capped = hplc.build_background(intensity, std, keep=0.25)[2]
        self.assertEqual(len(capped), 15, "60 candidates capped at a quarter")
        self.assertLess(len(capped), len(wide), "the ceiling must bite")
        self.assertTrue(set(capped) <= set(wide), "capping only removes the worst scores")
        # a ceiling of 100 % is no ceiling at all: the knee still has the last word
        numpy.testing.assert_array_equal(hplc.build_background(intensity, std, keep=1.0)[2], wide)

    def test_build_background_reports_the_three_metrics(self):
        """The two metrics that do not decide are returned beside the one that does.

        They are written next to the background so the cut can be audited: if the three
        rank the candidates alike, the frames retained owe nothing to the metric that
        happened to be chosen. Only their shape and alignment are checked here.
        """
        rng = numpy.random.default_rng(23)
        intensity = 100.0 * (1 + 0.002 * rng.standard_normal((80, 150)))
        std = 0.002 * intensity
        intensity[7] *= 1 + 0.02 * numpy.linspace(1, -1, 150)   # one curve bent out of shape
        bg = hplc.build_background(intensity, std)
        for name in ("zscore", "zscore_rms", "longest_streak"):
            self.assertEqual(getattr(bg.scores, name).shape, bg.candidates.shape,
                             f"{name} must carry one value per candidate")
        self.assertTrue(numpy.issubdtype(bg.scores.longest_streak.dtype, numpy.integer),
                        "a longest run is a count of bins")
        self.assertTrue((bg.scores.longest_streak <= intensity.shape[1]).all(),
                        "a run cannot be longer than the curve")
        numpy.testing.assert_allclose(bg.scores.zscore,
                                      hplc.background_zscore(intensity[bg.candidates],
                                                             *self._reference(intensity, std)),
                                      rtol=1e-12, err_msg="zscore must be the score that decides")
        self.assertTrue((bg.scores.zscore_rms >= bg.scores.zscore).all(),
                        "the L2 norm of a vector is never below its L1 mean")
        # the bent curve must stand out on every metric, whatever each one cuts at
        worst = numpy.flatnonzero(bg.candidates == 7)
        for name in ("zscore", "zscore_rms"):
            self.assertEqual(int(numpy.argmax(getattr(bg.scores, name))), int(worst[0]),
                             f"{name} must rank the bent curve worst")

    @staticmethod
    def _reference(intensity, std):
        "the SVD fundamental `build_background` scores against, and the matching σ"
        U, S, V = numpy.linalg.svd(intensity.T, full_matrices=False)
        return numpy.median(V[0]) * S[0] * U[:, 0], std

    def test_threshold_at_draws_the_cut(self):
        """The level written as a horizontal line puts exactly `nkeep` curves below it.

        It is what lets the three metrics of `5_background` be compared on one plot:
        each is cut where it would retain as many curves as the selection did.
        """
        values = numpy.array([5.0, 1.0, 4.0, 2.0, 3.0])
        for nkeep in range(1, 6):
            level = hplc.threshold_at(values, nkeep)
            self.assertEqual((values <= level).sum(), nkeep,
                             f"{nkeep} curves must sit at or below the line")
        self.assertTrue(numpy.isnan(hplc.threshold_at(values, 0)), "no curve, no line")
        self.assertTrue(numpy.isnan(hplc.threshold_at(values, 6)), "more curves than there are")

    def test_threshold_on_a_tied_metric_stops_separating(self):
        """Ties put more curves on the line than it is meant to retain.

        A count of bins — the CorMap statistic — takes few distinct values, so its line
        is at the right level yet no longer cuts cleanly. That staircase is the reason
        the score was replaced, and showing it is the point of plotting it.
        """
        values = numpy.array([1, 2, 2, 2, 2, 3], dtype=numpy.int32)
        level = hplc.threshold_at(values, 3)
        self.assertEqual(level, 2, "the third smallest value is still 2")
        self.assertEqual((values <= level).sum(), 5, "but five curves sit on or below it")

    def test_knee_is_not_invariant_to_squaring_the_score(self):
        """A monotone change of the score leaves the ranking alone but moves the cut.

        Squaring every score reorders nothing, so an L1 norm and an L2 one rank the
        candidates alike; the knee nevertheless answers differently, because it reads
        the *shape* of the sorted curve and not the order. Which way it moves is not
        predictable run by run, so only the non-invariance is pinned here — the
        direction, looser in the median, is a corpus measurement and lives in the
        docstring. This matters because an L2 score is the natural thing to reach for
        when the selection looks too permissive, and it is not a lever on that.
        """
        rng = numpy.random.default_rng(3)
        score = numpy.abs(rng.standard_normal((400, 200))).mean(axis=-1)
        plain = hplc.knee_threshold(score)
        squared = hplc.knee_threshold(score ** 2)
        numpy.testing.assert_array_equal(numpy.argsort(score), numpy.argsort(score ** 2),
                                         "squaring must not change the ranking")
        self.assertNotEqual(plain.sum(), squared.sum(), "yet it does change the cut")


class TestFractionWidth(unittest.TestCase):
    """Where `search_peaks` cuts the flanks of a peak, which PEAK_REL_HEIGHT sets."""

    @staticmethod
    def chromatogram(centre=100.0, width=10.0, nframes=200, npt=200, amplitude=50.0, seed=6):
        "A single gaussian elution on a flat solvent, quiet enough to be unambiguous"
        rng = numpy.random.default_rng(seed)
        q = numpy.linspace(0.05, 4.0, npt)
        frames = numpy.arange(nframes, dtype=numpy.float64)
        elution = numpy.exp(-0.5 * ((frames - centre) / width) ** 2)
        solute = numpy.exp(-(q / 0.3) ** 2)          # scatters at low q only
        intensity = numpy.ones((nframes, npt)) + amplitude * elution[:, None] * solute
        return q, intensity * (1 + 1e-4 * rng.standard_normal((nframes, npt)))

    def test_cuts_at_half_maximum(self):
        "The default width is the FWHM of the elution, 2.355 sigma"
        q, intensity = self.chromatogram(width=10.0)
        labels, count = hplc.search_peaks(intensity, q)
        self.assertEqual(count, 1, "one elution, one fraction")
        found = numpy.flatnonzero(labels == 1)
        self.assertAlmostEqual(found.size, 2.355 * 10.0, delta=3.0,
                               msg="the fraction is not the width at half maximum")

    def test_rel_height_widens_the_fraction(self):
        "Cutting lower on the flank takes more frames in, monotonically"
        q, intensity = self.chromatogram(width=10.0)
        sizes = [(hplc.search_peaks(intensity, q, rel_height=r)[0] == 1).sum()
                 for r in (0.5, 0.7, 0.9)]
        self.assertEqual(sizes, sorted(sizes), f"not monotonic: {sizes}")
        self.assertLess(sizes[0], sizes[-1], "rel_height has no effect")

    def test_default_is_the_constant(self):
        "The signature follows PEAK_REL_HEIGHT rather than a value of its own"
        q, intensity = self.chromatogram(width=10.0)
        numpy.testing.assert_array_equal(
            hplc.search_peaks(intensity, q)[0],
            hplc.search_peaks(intensity, q, rel_height=hplc.PEAK_REL_HEIGHT)[0])


def suite():
    loader = unittest.defaultTestLoader.loadTestsFromTestCase
    testsuite = unittest.TestSuite()
    testsuite.addTest(loader(TestAnalysis))
    testsuite.addTest(loader(TestBackgroundSelection))
    testsuite.addTest(loader(TestFractionWidth))
    return testsuite


if __name__ == "__main__":
    if "--regenerate" in sys.argv:
        generate_reference()
        print(f"Reference regenerated in {REFERENCE}")
    else:
        runner = unittest.TextTestRunner()
        if not runner.run(suite()).wasSuccessful():
            sys.exit(1)
