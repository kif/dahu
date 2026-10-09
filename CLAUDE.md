# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Plugins have their own guidance: `plugins/CLAUDE.md` (all plugins) and one file per beamline directory when it has specific constraints (e.g. `plugins/bm29/CLAUDE.md`).

## What dahu is

A lightweight, plugin-based data-analysis engine: a JSON-RPC server exposed over Tango (`dahu-server`), used at ESRF beamlines. A plugin takes one JSON-serializable dict as input and produces one as output; each execution is a `Job` running in its own thread.

## Build, install, test

Build system is **meson-python** (`pyproject.toml`, `meson.build`). Every Python file that has to be installed must be listed in the `meson.build` of its directory; a file missing there is absent from wheels and from the editable build (which fails if a listed file does not exist).

```bash
pip install .                                   # regular install
pip install --no-build-isolation -e .           # editable (needs meson-python, meson, ninja installed)

python run_tests.py                             # builds the project (bootstrap.py) then runs dahu.test.suite
python run_tests.py --installed                 # test the installed package instead of building
python -m dahu.test.test_all                    # core tests on an installed/editable dahu
python -m unittest dahu.test.test_job.TestJob.test_abort        # a single core test
```

- In editable mode the factory does not find `plugins/` (it looks next to `src/dahu`): set `DAHU_PLUGINS=$PWD/plugins` (`os.pathsep`-separated list) when running tests or `dahu-reprocess`.
- CI (`.github/workflows/python-package.yml`, Python 3.11–3.14): `flake8 . --select=E9,F63,F7,F82` must be clean (syntax errors / undefined names), then `python run_tests.py`.
- The core test suite only covers the kernel (`src/dahu`).

## Kernel architecture (`src/dahu`)

- `plugin.py` — `Plugin` base class: empty constructor, then `setup()` (reads/sanitizes `self.input`), `process()`, `teardown()` (fills `self.output`), optional `abort()`. `log_error(..., do_raise=True)` raises to fail the job; `wait_for(job_id)` synchronizes on another job. `plugin_from_function` turns a stateless function into a plugin.
- `factory.py` — `plugin_factory` (singleton `Factory`) maps a **lower-cased fully-qualified name** (`bm29.subtractbuffer`) to a class. On a cache miss it imports the module named by everything before the last dot from the plugin directories (`$DAHU_PLUGINS` first, then `dahu/plugins`); importing that module must call `register(...)`. `optional_plugin(fqn)` records why a plugin failed to load, reported instead of `None`.
- `job.py` — `Job(Thread)` instantiates the plugin, runs setup/process/teardown, stores status, and serializes input/output JSON to the working directory.
- `server.py` / `app/tango_server.py` — the Tango device; `app/reprocess.py` (`dahu-reprocess file.inp …`) replays saved job inputs offline and writes `NNNNN_<Plugin>.inp/.out` in a temporary work directory.
