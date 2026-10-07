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

## Pipeline

Chained jobs, each reading the HDF5/NeXus file written by the previous one:

1. `integrate.py` — `bm29.integratemultiframe`: azimuthal integration (pyFAI) of every frame of a sample-changer exposure → spottiness/meniscus detection → renormalization on the linear fit of the beam-stop diode (its variance is added) → CorMap selection of equivalent frames → average via pyFAI `IntegrateResult.union()`. The averaged curve and its accumulators (`sum_signal`, `sum_normalization`, `sum_normalization2`, `sum_variance`, `count`, `error_model` attribute) are stored in `4_time_average`. In `hplc_mode` it stops after `1_integration`.
2. `subtracte.py` — `bm29.subtractbuffer`: `NexusJuice.read()` loads integrated files; buffers equivalent by CorMap are averaged by rebuilding `Integrate1dResult`s from the accumulators (`NexusJuice.to_result()`) and merging them with `union()`; the subtraction is then done on 1D curves (errors summed quadratically) before the SAXS analysis.
3. `hplc.py` — `bm29.hplc`: rebuilds the chromatogram from frame-per-frame integrations, NMF, fractions, then the SAXS analysis per fraction.
4. `analysis.py` — SAXS analysis shared by subtract and HPLC: `saxs_analysis()` runs Guinier → dimensionless Kratky → invariants (Rambo-Tainer, Porod) → BIFT (freesas), each step as an `NXprocess` group, and stops after Guinier if no Guinier region is found.
5. `common.py` — shared types (`Sample`, `Ispyb`, `SequenceIndex`), compression settings `cmp_int`/`cmp_float` (hdf5plugin), integrator cache; `nexus.py` — NeXus writer.

`integrate.py` writes the file layout that `subtracte.NexusJuice.read()` and `hplc.NexusJuice.read()` parse: change them together.

## NeXus conventions

- Process groups are named `f"{index}_{name}"` with a `sequence_index` dataset taken from a `SequenceIndex` counter: keep names and indices consistent, no duplicated index.
- `NXdata` `signal`/`axes` attributes are relative dataset names, and every group sets its `default` attribute.

## External services

- ISPyB (`ispyb.py`) is contacted only when `ispyb.url` is set in the job input.
- iCAT (`icat.py`, `send_icat`) only proceeds when the gallery path contains a `PROCESSED_DATA` or `processed` directory.
- memcached (`memcached.py`) targets `localhost`.

To reprocess production job files offline without touching those services, remove `ispyb.url` (and `pyarch`) from the inputs and write outputs outside any `PROCESSED_DATA`/`processed` tree.
