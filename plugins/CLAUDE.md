# Plugins

Guidance common to all plugins. Beamline directories may add their own `CLAUDE.md` with stricter constraints.

- One module or package per beamline; they are independent from each other, do not import across beamlines. `example.py` shows class-based (`@register` on a `Plugin` subclass) and function-based plugins.
- A plugin is found by its lower-cased fully-qualified name `<module>.<plugin>`: the factory imports `<module>` (a `.py` file or a package `__init__.py` in this directory) and that import must `register` the plugin. In a package, wrap each registration in `with optional_plugin("<module>.<plugin>"):` so that a missing optional dependency disables only that plugin.
- Every installed file has to be listed in the `meson.build` of its directory (and a new sub-package needs a `subdir()` in the parent `meson.build`).
- Beamline-specific dependencies go in `<beamline>/requirements.txt`, not in the dahu core dependencies.
- Inputs and outputs are JSON: anything put in `self.output` must be serializable (numpy values are handled by `NumpyEncoder` when the job is saved). Keys of the input listed in `REPROCESS_IGNORE` are dropped by `dahu-reprocess`.
- Plugins are not covered by `run_tests.py`; only `bm29/test` exists so far.
