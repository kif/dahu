Installation of Dahu
====================

Dahu is a Python (-3) project, the only dependency is `numpy` but `Tango` is needed for the server part.

.. code-block:: shell

    python setup.py build test bdist_wheel
    pip install  dist/dahu*.whl


Dahu can also be installed as debian package with automatic packaging from the `build-deb.sh` tool.
Documentation can be generated with `python setup.py build_doc`.

Devlopment of Dahu
==================

`Dahu` is an Open source project under MIT license available at:
https://github.com/kif/dahu

There is a minimal test-suite which ensures the dynamic loading works as expected.
The kernel of dahu is tested, but plugins are not.


Development of plugins
======================

Plugins are regrouped into the `plugin` directory.
There is one directory per beamline and several plugins per beamline.
Each beamline is independant so no interference are expected.

Dahu can be tested on plugins in an alternative directory using the `$DAHU_PLUGINS` environment variable.

Different examples of plugins are provided in the plugins/example.py file, either based on a function or a class with the `@register` decorator.

