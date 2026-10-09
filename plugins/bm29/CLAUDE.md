# BM29 BioSAXS plugins

## Dependencies and tests

```bash
pip install -r plugins/bm29/requirements.txt     # numpy, scipy, pyFAI, h5py, hdf5plugin, freesas, ...
DAHU_PLUGINS=$PWD/plugins python -m dahu.plugins.bm29.test.test_analysis
python -m unittest dahu.plugins.bm29.test.test_analysis.TestAnalysis.test_subtract
python -m dahu.plugins.bm29.test.test_analysis --regenerate      # only when an output change is intended
```

`hplc.py` needs a freesas recent enough to provide `freesas.containers.UVJuice` and `freesas.plot.hplc_plot`.

`test/test_analysis.py` runs HPLC and SubtractBuffer on synthetic photon-counting data (sphere of radius 3 nm, Poisson noise, file reading mocked) and compares a summary of the produced HDF5 to `test/reference_analysis.json`, ignoring dates, versions and environment-dependent fields. After `--regenerate`, review the diff of the JSON: it must only contain the intended changes.

The reference also drifts on its own when freesas is upgraded: the BIFT of the `02s-12s`
fraction fits badly on purpose, so its Powell descent is version-sensitive. Expect that
fraction, and no other, to move. It mocks `NexusJuice`, so neither the diode, the
concatenation nor the peak search are covered: those are exercised against real files,
not by this suite.

## Pipeline

Chained jobs, each reading the HDF5/NeXus file written by the previous one:

1. `integrate.py` — `bm29.integratemultiframe`: azimuthal integration (pyFAI) of every frame → spottiness/meniscus detection → in sample-changer mode only, renormalization on the linear fit of the beam-stop diode (its variance is added), CorMap selection of equivalent frames and average via `IntegrateResult.union()` in `4_time_average`. In `hplc_mode` it stops after `1_integration`. Accumulators are written per frame in `1_integration/accumulators` and, in SC mode, again in `2_renormalize/accumulators`: `sum_signal`, `sum_normalization`, `sum_normalization2`, `sum_variance_azimuthal`, `sum_variance_poisson`, `count`, all ZFP-compressed.
2. `subtracte.py` — `bm29.subtractbuffer`: `NexusJuice.read()` loads integrated files; buffers equivalent by CorMap are averaged by rebuilding `Integrate1dResult`s from the accumulators (`NexusJuice.to_result()`) and merging them with `union()`; the subtraction is then done on 1D curves (errors summed quadratically) before the SAXS analysis.
3. `hplc.py` — `bm29.hplc`: see below.
4. `analysis.py` — SAXS analysis shared by subtract and HPLC: `saxs_analysis()` runs Guinier → dimensionless Kratky → invariants (Rambo-Tainer, Porod) → BIFT (freesas), each step as an `NXprocess` group. It stops right after Guinier when no region is found **or when the retained fit starts past `GUINIER_QRG_LIMIT`**, in which case `2_Guinier_analysis/invalid` says why.
5. `common.py` — shared types (`Sample`, `Ispyb`, `SequenceIndex`), compression settings `cmp_int`/`cmp_float` (hdf5plugin), integrator cache; `nexus.py` — NeXus writer.

`integrate.py` writes the file layout that `subtracte.NexusJuice.read()` and `hplc.NexusJuice.read()` parse: change them together.

### The HPLC plugin, step by step

`hplc.NexusJuice.read()` collects everything `integrate.py` leaves behind, including the
accumulators, the diode, the ring current, the spottiness and the LImA file the frames
are an external link to. `NexusJuice.concatenate()` then merges the juices of all the
input files **in frame order**: nothing guarantees the jobs finish in order, and the
frame index is the only reliable ordering. Series are scattered into arrays of
`max(frame) + 1` entries, so a frame no file provides stays at zero.

The NeXus output is then:

```
0_measurement    external links to the integrate files
1_renormalize    diode/ (raw, smooth, smooth_errors, frame_idx, timestamps)
                 result/ (I, errors, q)        accumulators/ (the six, corrected)
2_chromatogram   SAXS/ (sum_q_range, sum, diode)   UV-Vis/   result/ (SAXS + UV, 0-1)
3_SVD            4_NMF            5_background            6_SEC_fractions
7_ISPyB
```

**1_renormalize.** The diode is smoothed and the curves are divided by the smoothed
value rather than by the raw reading. Measured over 177 runs of a season: the noise of
the diode goes into the curves one for one (correlation 0.86, slope 1.07), and taking it
out removes some 75 % of the scatter of the chromatogram. The uncertainty of the
smoothed diode is injected into both variances as `sum_signal²·var_d/d²`, exactly as
`integrate.py` does in the sample-changer path. The accumulators are corrected for the
renormalization — `sum_normalization` scaled by the inverse factor, `sum_normalization2`
by its square — and are what everything downstream merges.

Mind that `sum_normalization2` only keeps its statistical meaning until a
renormalization happens; past that it is propagated for the sake of pyFAI's machinery
alone. Downstream, use `sem` and never `std`.

**2_chromatogram.** `SAXS` carries the chromatogram summed over `CHROMATOGRAM_QRANGE` as
its signal, with the whole-range sum and the diode as auxiliary signals. `result`
overlays it with the UV-Vis channels, all rescaled between 0 and 1. The UV counters are
read from the BLISS master file, found by walking up from the LImA file the frames link
to (`find_bliss_master`): that file does not exist yet when `integrate.py` runs, but it
does by now, and its counters share the clock of the frames. The `.dat` written by the
spectrometer is the fallback, and it needs `uv_offset` because its clock is its own.
Every access to the BLISS master is guarded and done without locking: the acquisition
may still hold it open for writing.

**5_background — selecting the buffer frames.** Three criteria, answering three
different questions, and none replaces another:

- `stationary_frames()` — *is this the solvent of the run?* The level over
  `BACKGROUND_QRANGE`, where the solute scatters almost nothing, identifies the liquid in
  the capillary. A frame off by more than `nsigma` robust standard deviations saw a
  bubble, another solvent or a glitch. Blind to the solute.
- `solute_frames()` — *is anything eluting?* The ratio of `SOLUTE_QRANGE` to
  `BACKGROUND_QRANGE` rises with the solute and does not care about the overall scale. It
  covers the whole elution, shoulders included, which matters: leaving the shoulders of a
  peak in a background is worse than leaving the whole peak, because they drag the SVD
  fundamental along without being obvious enough for cormap to reject them.
- CorMap against the fundamental of the SVD, on what the two gates left — the historical
  criterion, kept for what they do not see: crystallites, parasitic scattering.

`build_background()` then averages `background_keep` of the survivors **through the
accumulators**, not as a mean of ratios: frames are weighted by their own normalization,
which matters as soon as the flux varies or the detector marks part of an image invalid
mid-run.

On the dataset the gates were built from, they took the bias of the background from
−0.76 % to +0.03 %, the high-q residual of the subtracted peak from +1.07 % to +0.37 %,
and the Porod volume from 49 to 91 nm³. Watch the fragmentation of `kept`: a buffer
selection scattered over dozens of two-frame blocks means the criteria no longer
discriminate, and is worth a look before trusting the result.

**6_SEC_fractions — finding the peaks.** `search_peaks()` runs `scipy.signal.find_peaks`
on the chromatogram summed over `CHROMATOGRAM_QRANGE`, with a prominence of
`PEAK_PROMINENCE` times the noise measured on that very chromatogram by
`estimate_noise()`. Regions are cut at the lowest point between two summits so that they
stay disjoint; taking the width at the base of each peak would nest a shoulder inside its
neighbour. Two safety nets: an elution still rising at the last frame has no summit and
is picked up from `solute_frames()` instead, and a flag confined to the first or last
`PEAK_EDGE_MARGIN` frames, or sitting on frames `stationary_frames()` rejects, is
discarded — those are start-up artefacts and solvent changes, not elutions.

Benchmarked over 174 runs against a physical ground truth (low-q excess above 8 robust
deviations): 4.46 regions per run, recall 0.95, precision 0.33 for the wavelet transform
it replaces, against 1.55, 0.98 and 0.82. The wavelet version was dominated on **both**
axes, so tolerating its false positives bought no recall.

Each fraction is merged through the accumulators as well, and the background is
subtracted **once at the end** rather than frame by frame: it is the same curve removed
from every frame, so its uncertainty is fully correlated and does not average down.
Subtracting first understates it by √n — 11 % on a fraction of 129 frames.

The zip beside the HDF5 holds ASCII `.dat` files: `buffer.dat` and
`fraction_<first>-<last>.dat` at the root, every frame under `frames/`. Everything there
is named by frame index; time stamps live in the HDF5, where they serve the comparison
with the UV-Vis chromatogram.

### Constants to revisit

All in `hplc.py` unless noted, each with the measurement behind it in its docstring.
They were tuned on BM29 data of 2026 and are the first things to look at if the beamline,
the detector or the elution conditions change.

| Constant | Value | What it is |
|---|---|---|
| `CHROMATOGRAM_QRANGE` | `(0.1, 1.0)` nm⁻¹ | Summed to build the chromatogram and to look for peaks. The bins above carry solvent and noise only; dropping them buys a factor 3 on the signal to noise. |
| `BACKGROUND_QRANGE` | `(1.5, 4.0)` nm⁻¹ | Where the solvent scatters alone: identifies the liquid. |
| `SOLUTE_QRANGE` | `(0.1, 0.5)` nm⁻¹ | Where the solute shows: its ratio to the above detects an elution. |
| `DIODE_FILTER_SIZE` | `11` frames | Width of the diode filter, overridden by the `diode_medfilt` input. Wider is not better: the diode drifts and steps, and at 31 frames the filter biased one run out of six by more than 3σ. |
| `DIODE_NOISE_LIMIT` | `1.0` % | Above it, a `log_warning` fires: normalizing on such a diode is meaningless. Four runs out of 177 hit it, all beam-off test scans. |
| `PEAK_PROMINENCE` | `10.0` × noise | How far a peak must stand out. Insensitive between 5 and 25: precision moves from 0.81 to 0.85, recall from 0.97 to 0.94. |
| `PEAK_EDGE_MARGIN` | `10` frames | A lone excess confined to either end is a start-up artefact. A truncated elution spans hundreds of frames and is kept. |
| `ZIP_COMPRESSION` | `ZIP_DEFLATED` | bzip2 and lzma give 5.9 and 4.4 MB against 7.0, but Windows does not read them and these files go to users. |
| `GUINIER_QRG_LIMIT` | `1.0` (`analysis.py`) | Largest q·Rg at which a fit may start. Past it Guinier and BIFT disagree by 52 % on Rg, against 7 % below 0.6. |

Two more are not constants but hardcoded at their call site, and deserve to become
inputs the day they need tuning: `nsigma = 5.0`, the width of both gates, and
`minimum_size = 10`, the narrowest peak worth a fraction.

Inputs of `bm29.hplc`, beyond the files: `nmf_components`, `diode_medfilt`,
`diode_filter` (one of `SMOOTHING_ALGORITHMS`), `uv_datafile`, `uv_offset`,
`background_keep`.

## NeXus conventions

- Process groups are named `f"{index}_{name}"`: keep the numbering contiguous and never duplicate an index. `1_renormalize` is written even when nothing is smoothed, precisely so that the numbering does not depend on the options.
- The `sequence_index` dataset comes from a separate `SequenceIndex` counter which sub-groups consume as well, so it drifts from the number in the name — in a typical HPLC run `3_SVD` carries index 4 and `7_ISPyB` index 28. Do not read one from the other.
- `NXdata` `signal`/`axes` attributes are relative dataset names, and every group sets its `default` attribute.

## External services

- ISPyB (`ispyb.py`) is contacted only when `ispyb.url` is set in the job input.
- iCAT (`icat.py`, `send_icat`) only proceeds when the gallery path contains a `PROCESSED_DATA` or `processed` directory.
- memcached (`memcached.py`) targets `localhost`.

To reprocess production job files offline without touching those services, remove `ispyb.url` (and `pyarch`) from the inputs and write outputs outside any `PROCESSED_DATA`/`processed` tree.
