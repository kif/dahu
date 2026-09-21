#!/usr/bin/env python3
#

"""Test suite for the plugin factory, mainly the isolation of the registration."""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "14/09/2026"
__status__ = "production"

import os
import shutil
import sys
import tempfile
import unittest

from ..factory import Factory, optional_plugin, plugin_factory
from . import utilstest

logger = utilstest.getLogger(__name__)

PACKAGE = "dahutestbl"  # name of the fake beamline package built on the fly

BROKEN_PLUGIN = '''"Plugin relying on a third party module which is not installed"
import module_which_does_not_exist  # noqa
from dahu.plugin import Plugin


class Broken(Plugin):
    "Plugin which can never be loaded"
'''

GOOD_PLUGIN = '''"Plugin without any exotic dependency"
from dahu.plugin import Plugin


class Good(Plugin):
    "Plugin which works"

    def process(self):
        self.output["result"] = 42
'''

INIT_PLUGIN = '''"Fake beamline used to validate the isolation of the registration"
from dahu.factory import optional_plugin, register

# The broken plugin comes first on purpose: it used to abort the import of the
# whole package, making all the other plugins of the beamline unavailable.
with optional_plugin("%(pkg)s.broken"):
    from .broken import Broken
    register(Broken, fqn="%(pkg)s.broken")

with optional_plugin("%(pkg)s.good"):
    from .good import Good
    register(Good, fqn="%(pkg)s.good")
''' % {"pkg": PACKAGE}


class TestFactory(unittest.TestCase):
    """Build a fake beamline package in a temporary directory and load it."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="dahu_test_factory_")
        pkgdir = os.path.join(cls.tmpdir, PACKAGE)
        os.mkdir(pkgdir)
        for filename, content in (("__init__.py", INIT_PLUGIN),
                                  ("broken.py", BROKEN_PLUGIN),
                                  ("good.py", GOOD_PLUGIN)):
            with open(os.path.join(pkgdir, filename), "w") as fd:
                fd.write(content)
        plugin_factory.add_directory(cls.tmpdir)

    @classmethod
    def tearDownClass(cls):
        "Remove every trace of the fake package from the (class-wide) factory"
        Factory.plugin_dirs.pop(os.path.abspath(cls.tmpdir), None)
        Factory.modules.pop(PACKAGE, None)
        sys.modules.pop(PACKAGE, None)
        for fqn in (f"{PACKAGE}.good", f"{PACKAGE}.broken"):
            Factory.registry.pop(fqn, None)
            Factory.unavailable.pop(fqn, None)
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_working_plugin_survives_broken_sibling(self):
        "A plugin must stay available when another one of the package fails to load"
        plugin = plugin_factory(f"{PACKAGE}.good")
        self.assertIsNotNone(plugin, "plugin registered despite its broken sibling")
        plugin.process()
        self.assertEqual(plugin.output["result"], 42, "plugin is usable")

    def test_broken_plugin_is_disabled(self):
        "A plugin which fails to load is disabled and the reason is recorded"
        self.assertIsNone(plugin_factory(f"{PACKAGE}.broken"), "plugin not registered")
        reason = Factory.unavailable.get(f"{PACKAGE}.broken")
        self.assertIsNotNone(reason, "reason of the failure is recorded")
        self.assertIn("module_which_does_not_exist", reason, "reason is explicit")

    def test_optional_plugin(self):
        "The context manager swallows the error, records it, then forgets it"
        fqn = "nowhere.nothing"
        try:
            with optional_plugin(fqn):
                raise RuntimeError("some failure")
            self.assertEqual(Factory.unavailable.get(fqn),
                             "RuntimeError: some failure", "failure recorded")
            with optional_plugin(fqn):
                pass
            self.assertNotIn(fqn, Factory.unavailable, "stale failure discarded")
        finally:
            Factory.unavailable.pop(fqn, None)


def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(TestFactory("test_working_plugin_survives_broken_sibling"))
    testSuite.addTest(TestFactory("test_broken_plugin_is_disabled"))
    testSuite.addTest(TestFactory("test_optional_plugin"))
    return testSuite


if __name__ == '__main__':
    mysuite = suite()
    runner = unittest.TextTestRunner()
    runner.run(mysuite)
