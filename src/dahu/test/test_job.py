#!/usr/bin/env python3
#

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "11/03/2026"
__status__ = "production"

import os
import unittest

from .. import job
from . import utilstest

logger = utilstest.getLogger(__name__)


class TestJob(unittest.TestCase):
    def test_plugin_from_function(self):
        self.called = False
        dico = {"plugin_name": "example.square",
                "x": 5}
        j = job.Job("example.square", dico)
        j.connect_callback(self.callback)
        logger.info(j)
        j.start()
        j.join()
        logger.info(j)
        if "error" in j.output_data:
            logger.error(os.linesep.join(j.output_data["error"]))

        logger.info(j.input_data)
        self.assertTrue(self.called)
        print(j.output_data)
        self.assertEqual(j.output_data["result"], 25, "result OK")

    def callback(self, *args, **kwargs):
        logger.info(f"callback actually called with  arguments {args} and kwargs {kwargs}")
        assert len(args) == 1
        self.called = True

    def test_callbacks_on_missing_plugin(self):
        "A job which fails to instanciate its plugin must still run its callbacks"
        self.called = False
        j = job.Job("nosuch.plugin", {})
        j.connect_callback(self.callback)
        j.start()  # synchronous: the thread is never started in this case
        self.assertEqual(j.status, j.STATE_FAILURE, "job ended in failure")
        self.assertTrue(self.called, "callback called despite the missing plugin")

    def test_clean_job_from_id(self):
        "Cleaning a job frees the plugin and leaves the data readable from disk"
        j = job.Job("example.square", {"x": 5})
        j.start()
        j.join()
        self.assertEqual(j.status, j.STATE_SUCCESS, "job succeeded")

        msg = job.Job.clean_job_from_id(j.id)
        self.assertEqual(msg, f"Job {j.id} cleaned", "job cleaned")
        self.assertTrue(j.data_on_disk, "input and output serialized on disk")
        for ext in (".inp", ".out"):
            self.assertTrue(os.path.isfile(j.data_on_disk + ext), f"{ext} file written")
        output = job.Job.getDataOutputFromId(j.id)
        self.assertEqual(output["result"], 25, "output still readable once cleaned")

        unknown = job.Job.clean_job_from_id(j.id + 1000)
        self.assertIn("Unable to retrieve", unknown, "unknown job does not raise")

    def test_clean_job_from_id_aliases(self):
        "The camelCase spellings are kept for backward compatibility"
        aliases = ("cleanJobfromId", "cleanJobfromID",
                   "cleanJobFromId", "cleanJobFromID")
        for name in aliases:
            self.assertTrue(hasattr(job.Job, name), f"Job.{name} exists")
            self.assertIs(getattr(job.Job, name).__func__,
                          job.Job.clean_job_from_id.__func__,
                          f"Job.{name} is clean_job_from_id")


def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(TestJob("test_plugin_from_function"))
    testSuite.addTest(TestJob("test_callbacks_on_missing_plugin"))
    testSuite.addTest(TestJob("test_clean_job_from_id"))
    testSuite.addTest(TestJob("test_clean_job_from_id_aliases"))
    return testSuite

if __name__ == '__main__':
    mysuite = suite()
    runner = unittest.TextTestRunner()
    runner.run(mysuite)
