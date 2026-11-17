#!/usr/bin/env python3
#

"""Test suite for the Tango device server.

The device server is tested without a Tango installation: PyTango is replaced by
a stand-in when it is missing, which is enough to exercise the plugin reporting.
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "14/09/2026"
__status__ = "production"

import inspect
import os
import sys
import types
import unittest

from ..factory import Factory
from . import utilstest

logger = utilstest.getLogger(__name__)

FAKE_FQN = "nowhere.nothing"
FAKE_REASON = "ModuleNotFoundError: No module named 'nothing'"


def stub_pytango():
    """Build a minimal stand-in for PyTango

    Only what is needed to import `dahu.server`: the two base classes and the
    constants used in the class bodies.

    :return: a module object to be put in `sys.modules`
    """
    tango = types.ModuleType("PyTango")

    class _Base:

        def __init__(self, *args, **kwargs):
            pass

    tango.LatestDeviceImpl = _Base
    tango.DeviceClass = _Base
    for name in ("DevString", "DevLong", "DevBoolean", "DevVoid",
                 "DevVarStringArray", "SCALAR", "READ", "READ_WRITE",
                 "DevState", "Util", "Database"):
        setattr(tango, name, type(name, (), {}))
    tango.DevFailed = type("DevFailed", (Exception,), {})
    return tango


class FakeAttribute:
    """Stand-in for a Tango attribute, enough for the `serialize` accessors"""

    def __init__(self, write_value=None):
        self.write_value = write_value
        self.value = None

    def get_write_value(self):
        return self.write_value

    def set_value(self, value):
        self.value = value


class TestServer(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            import PyTango  # noqa: F401
        except ImportError:
            logger.info("PyTango is missing, using a stand-in for it")
            sys.modules["PyTango"] = stub_pytango()
            cls.stubbed = True
        else:
            cls.stubbed = False
        from ..server import DahuDS
        cls.DahuDS = DahuDS
        # the methods under test only use get_name(), no need for a real device
        cls.device = types.SimpleNamespace(get_name=lambda: "DahuDS")

    @classmethod
    def tearDownClass(cls):
        if cls.stubbed:
            sys.modules.pop("dahu.app.tango_server", None)
            sys.modules.pop("dahu.server", None)
            sys.modules.pop("PyTango", None)

    def test_summarize(self):
        "A plugin without docstring must not break the listing"
        from ..server import summarize
        self.assertEqual(summarize(None), "no documentation", "no docstring at all")
        self.assertEqual(summarize("   "), "no documentation", "blank docstring")
        self.assertEqual(summarize("\n\n  Some plugin\n  more text\n"),
                         "Some plugin", "first non empty line")

    def test_list_plugins(self):
        "Plugins which failed to load are listed with the reason"
        listing = self.DahuDS.listPlugins(self.device)
        self.assertNotIn(FAKE_FQN, listing, "nothing reported when all is fine")
        Factory.unavailable[FAKE_FQN] = FAKE_REASON
        try:
            listing = self.DahuDS.listPlugins(self.device)
        finally:
            Factory.unavailable.pop(FAKE_FQN, None)
        self.assertIn(FAKE_FQN, listing, "disabled plugin is listed")
        self.assertIn(FAKE_REASON, listing, "reason is displayed")

    def test_init_plugin(self):
        "initPlugin tells why a plugin is missing instead of just `None`"
        Factory.unavailable[FAKE_FQN] = FAKE_REASON
        try:
            msg = self.DahuDS.initPlugin(self.device, FAKE_FQN)
        finally:
            Factory.unavailable.pop(FAKE_FQN, None)
        self.assertIn(FAKE_REASON, msg, "reason of the failure is reported")

        msg = self.DahuDS.initPlugin(self.device, "nosuch.plugin")
        self.assertIn("no such plugin", msg, "unknown plugin reported as such")

        msg = self.DahuDS.initPlugin(self.device, "example.cube")
        self.assertTrue(msg.startswith("Plugin loaded"), "valid plugin still loads")

    def test_abort(self):
        "The abort command returns a boolean, as declared in cmd_list"
        result = self.DahuDS.abort(self.device, 999999)
        self.assertIsInstance(result, bool, "a DevBoolean is returned")
        self.assertFalse(result, "nothing to abort")

    def test_no_dead_command(self):
        "Every public method of the device must be reachable from Tango"
        from ..server import DahuDSClass
        internal = {"init_device", "delete_device", "always_executed_hook",
                    "read_attr_hardware", "get_name", "process_job",
                    "process_event", "finished_processing", "statistics"}
        for name, method in inspect.getmembers(self.DahuDS, inspect.isfunction):
            if name.startswith(("_", "read_", "write_")) or name in internal:
                continue
            if method.__module__ != "dahu.server":
                continue
            self.assertIn(name, DahuDSClass.cmd_list, f"{name} is declared as a command")

    def test_command_line(self):
        "The options of Tango must reach it untouched, not be eaten by argparse"
        from ..app.tango_server import parse
        argv = ["dahu_server", "dahu", "-ORBendPoint", "giop:tcp::10001",
                "-nodb", "-dlist", "id00/dahu/1", "-v4"]
        options, tangoParam = parse(argv)
        self.assertEqual(tangoParam, ["DahuDS"] + argv[1:], "options passed through")
        self.assertFalse(options.debug, "-dlist did not switch the debug mode on")
        self.assertNotIn("ist", os.path.basename(options.dahu_log),
                         "-dlist did not hijack the log directory")

        options, tangoParam = parse(["dahu_server", "-d", "-l", "/tmp/log", "dahu", "-nodb"])
        self.assertTrue(options.debug, "the options of dahu are still parsed")
        self.assertEqual(options.dahu_log, "/tmp/log", "log directory honoured")
        self.assertEqual(tangoParam, ["DahuDS", "dahu", "-nodb"], "instance then options")

    def test_serialize_attribute(self):
        "The serialize attribute must honour the value which is written"
        device = types.SimpleNamespace(_serialize=None)
        for written in (True, False):
            attr = FakeAttribute(written)
            self.DahuDS.write_serialize(device, attr)
            self.assertIs(device._serialize, written, f"serialize set to {written}")
            self.DahuDS.read_serialize(device, attr)
            self.assertIs(attr.value, written, f"serialize read back as {written}")


def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(TestServer("test_summarize"))
    testSuite.addTest(TestServer("test_list_plugins"))
    testSuite.addTest(TestServer("test_init_plugin"))
    testSuite.addTest(TestServer("test_abort"))
    testSuite.addTest(TestServer("test_no_dead_command"))
    testSuite.addTest(TestServer("test_command_line"))
    testSuite.addTest(TestServer("test_serialize_attribute"))
    return testSuite


if __name__ == '__main__':
    mysuite = suite()
    runner = unittest.TextTestRunner()
    runner.run(mysuite)
