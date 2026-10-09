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


def suite():
    loader = unittest.defaultTestLoader.loadTestsFromTestCase
    testsuite = unittest.TestSuite()
    testsuite.addTest(loader(TestAnalysis))
    return testsuite


if __name__ == "__main__":
    if "--regenerate" in sys.argv:
        generate_reference()
        print(f"Reference regenerated in {REFERENCE}")
    else:
        runner = unittest.TextTestRunner()
        if not runner.run(suite()).wasSuccessful():
            sys.exit(1)
