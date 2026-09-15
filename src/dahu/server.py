#!/usr/bin/env python3

"""
Data Analysis RPC server over Tango:

Tango device server
"""
__author__ = "Jérôme Kieffer"
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "11/03/2026"
__status__ = "production"
__docformat__ = 'restructuredtext'

import logging
import os
import sys
import threading
import time
from queue import Queue

import PyTango

from .job import Job, plugin_factory

logger = logging.getLogger("dahu.server")
# set loglevel at least at INFO
if logger.getEffectiveLevel() > logging.INFO:
    logger.setLevel(logging.INFO)

try:
    from rfoo.utils import rconsole
    rconsole.spawn_server()
except ImportError:
    logger.debug("No socket opened for debugging -> please install rfoo")


def summarize(docstring):
    """Extract a one line description out of a docstring

    :param docstring: the docstring of a plugin, possibly None
    :return: the first non empty line of it
    """
    for line in (docstring or "").splitlines():
        line = line.strip()
        if line:
            return line
    return "no documentation"


class DahuDS(PyTango.LatestDeviceImpl):
    """
    Tango device server launcher for Dahu server.
    """
    def __init__(self, cl, name):
        PyTango.LatestDeviceImpl.__init__(self, cl, name)
        self.init_device()
        self.job_queue = Queue()  # queue containing jobs to process
        self.event_queue = Queue()  # queue containing finished jobs
        self.processing_lock = threading.Semaphore()
        self.stat_lock = threading.Semaphore()
        self.last_stats = "No statistics collected yet, please use the 'collectStatistics' method first"
        self.last_failure = -1
        self.last_success = -1
        self.statistics_threads = None
        self._serialize = False
#         self._ncpu_sem = threading.Semaphore(multiprocessing.cpu_count())
        # start the two threads related to queues: process_job and event_queue
        t2 = threading.Thread(target=self.process_job)
        t2.start()
        t1 = threading.Thread(target=self.process_event)
        t1.start()


    def get_name(self):
        """Returns the name of the class"""
        return self.__class__.__name__

    def delete_device(self):

        logger.debug(f"[Device delete_device method] for device {self.get_name()}")

    def init_device(self):
        logger.debug(f"In {self.get_name()}.init_device()")

        self.set_state(PyTango.DevState.ON)
        self.get_device_properties(self.get_device_class())
        self.set_change_event("jobSuccess", True, False)
        self.set_change_event("jobFailure", True, False)
        self.set_change_event("statisticsCollected", True, False)

    def always_executed_hook(self):
        pass

    def read_attr_hardware(self, data):
        logger.debug(f"In {self.get_name()}.read_attr_hardware()")

    def read_jobSuccess(self, attr):
        attr.set_value(self.last_success)

    def read_jobFailure(self, attr):
        attr.set_value(self.last_failure)

    def read_statisticsCollected(self, attr):
        attr.set_value(self.last_stats)

    def read_serialize(self, attr):
        attr.set_value(bool(self._serialize))

    def write_serialize(self, attr):
        self._serialize = bool(attr.get_write_value())

    def getJobState(self, jobId):
        return Job.getStatusFromID(jobId)

    def cleanJob(self, jobId):
        return Job.clean_job_from_id(jobId)

    def listPlugins(self):
        """
        List all plugin currently loaded .... with a brief description
        """
        logger.debug(f"In {self.get_name()}.listPlugins")
        res = ["List of all plugin currently loaded (use initPlugin to loaded additional plugins):"]
        for name in sorted(plugin_factory.registry):
            res.append(f" {name} : {summarize(plugin_factory.registry[name].__doc__)}")
        unavailable = plugin_factory.unavailable
        if unavailable:
            res.append("")
            res.append("Plugins which failed to be loaded:")
            for name in sorted(unavailable):
                res.append(f" {name} : {unavailable[name]}")
        return os.linesep.join(res)

    def initPlugin(self, name):
        """
        Creates a job with the given plugin
        """
        logger.debug(f"In {self.get_name()}.initPlugin({name})")
        plugin = err = None
        try:
            plugin = plugin_factory(name)
        except Exception as error:
            err = f"{type(error).__name__}: {error}"
            logger.error(f"plugin {name} failed to be instanciated, {err}")
        if plugin is None:
            if err is None:
                # the plugin may have been disabled at registration time
                err = plugin_factory.unavailable.get(name.lower(), "no such plugin")
            return f"Plugin not found: {name}, {err}"
        else:
            return f"Plugin loaded: {name}{os.linesep}{plugin.__doc__}"

    def abort(self, jobId):
        """
        Aborts a job

        The plugin has to honour it: it is expected to check `is_aborted`.

        @param  jobId: ID of the job to stop
        @return: True if the job was running and has been asked to stop
        """
        logger.debug(f"In {self.get_name()}.abort({jobId})")
        return Job.abort_job_from_id(jobId)

    def startJob(self, argin):
        """
        Starts a job

        @param argin: 2-list [<Dahu plugin to execute>, <JSON serialized dict>]
        @return: jobID which is an int (-1 for error)
        """
        logger.debug(f"In {self.get_name()}.startJob()")
        name, data_input = argin[:2]
        if data_input.strip() == "":
            return -1
        job = Job(name, data_input)
        if job is None:
            return -1
        self.job_queue.put(job)
        return job.id

    def process_job(self):
        """
        Process all jobs in the queue.
        """
        while True:
            job = self.job_queue.get()
            job.connect_callback(self.finished_processing)
            job.start()
            if self._serialize:
                job.join()

    def finished_processing(self, job):
        """
        callback: when processing is done

        @param job: instance of dahu.job.Job
        """
        logger.debug(f"In {self.get_name()}.finished_processing id:{job.id} ({job.status})")
#         self._ncpu_sem.release()
        job.clean(wait=False)
        if job.status == job.STATE_SUCCESS:
            self.last_success = job.id
        else:
            sys.stdout.flush()
            sys.stderr.flush()
            self.last_failure = job.id
        self.job_queue.task_done()
        self.event_queue.put(job)

    def process_event(self):
        """
        process finished jobs on the tango side (issue with tango locks)
        """

        while True:
            job = self.event_queue.get()
            if job.status == job.STATE_SUCCESS:
                self.push_change_event("jobSuccess", job.id)
            else:
                self.push_change_event("jobFailure", job.id)

# TODO one day
#    def getRunning(self):
#        """
#        retrieve the list of plugins currently under execution (with their plugin-Id)
#        """
#        return EDStatus.getRunning()
#
#    def getSuccess(self):
#        """
#        retrieve the list of plugins finished with success (with their plugin-Id)
#        """
#        return EDStatus.getSuccess()
#
#    def getFailure(self):
#        """
#        retrieve the list of plugins finished with failure (with their plugin-Id)
#        """
#        return EDStatus.getFailure()

    def collectStatistics(self):
        """
        Retrieve some statistics on all Dahu-Jobs
        @return: a page of information about Dahu-jobs
        """
        self.statistics_threads = threading.Thread(target=self.statistics)
        self.statistics_threads.start()

    def statistics(self):
        """
        retrieve some statistics about past jobs.
        """
        with self.stat_lock:
            fStartStat = time.time()
            self.last_stats = Job.stats()
            self.last_stats += os.linesep + f"Statistics collected on {time.asctime()}, the collect took: {time.time() - fStartStat:.3f}s"
            self.push_change_event("statisticsCollected", self.last_stats)

    def getStatistics(self):
        """
        just return statistics previously calculated
        """
        if self.statistics_threads:
            self.statistics_threads.join()
        return self.last_stats

    def getJobOutput(self, jobId):
        """
        Retrieve XML output form a job
        @param jobId: name of the job
        @return: output from a job
        """
        return Job.getDataOutputFromId(jobId, as_JSON=True)

    def getJobInput(self, jobId):
        """
        Retrieve input from a job as JSON string
        @param jobId: identifier of the job (int)
        @return: JSON serialized input from a job
        """
        return Job.getDataInputFromId(jobId, as_JSON=True)

    def getJobError(self, jobId):
        """
        Retrieve error message from a job as a string
        @param jobId: identifier of the job (int)
        @return: Error message
        """
        return Job.getErrorFromId(jobId)

    def waitJob(self, jobId):
        """
        Wait for a job to be finished and returns the status.
        May cause Tango timeout if too slow to finish ....
        May do polling to wait the job actually started

        @param jobId: identifier of the job (int)
        @return: status of the job
        """
        res = Job.synchronize_job(jobId)
        i = 0
        while res == Job.STATE_UNINITIALIZED:
            if i > 10:
                break
            i += 1
            time.sleep(0.1)
            res = Job.synchronize_job(jobId)
        return res


class DahuDSClass(PyTango.DeviceClass):
    #    Class Properties
    class_property_list = {
        }

    #    Device Properties
    device_property_list = {
        'plugins_directory':
            [PyTango.DevString,
            "Dahu plugins directory",
            [] ],
        }


    #    Command definitions
    cmd_list = {
        'startJob': [[PyTango.DevVarStringArray, "[<Dahu plugin to execute>, <JSON serialized dict>]"], [PyTango.DevLong, "job id"]],
        'abort': [[PyTango.DevLong, "job id"], [PyTango.DevBoolean, ""]],
        'getJobState': [[PyTango.DevLong, "job id"], [PyTango.DevString, "job state"]],
        'initPlugin': [[PyTango.DevString, "plugin name"], [PyTango.DevString, "Message"]],
        'cleanJob':[[PyTango.DevLong, "job id"], [PyTango.DevString, "Message"]],
        'collectStatistics':[[PyTango.DevVoid, "nothing needed"], [PyTango.DevVoid, "Collect some statistics about jobs within Dahu"]],
        'getStatistics':[[PyTango.DevVoid, "nothing needed"], [PyTango.DevString, "Retrieve statistics about Dahu-jobs"]],
        'getJobOutput': [[PyTango.DevLong, "job id"], [PyTango.DevString, "<JSON serialized dict>"]],
        'getJobInput': [[PyTango.DevLong, "job id"], [PyTango.DevString, "<JSON serialized dict>"]],
        'getJobError': [[PyTango.DevLong, "job id"], [PyTango.DevString, "Error message"]],
        'listPlugins': [[PyTango.DevVoid, "nothing needed"], [PyTango.DevString, "prints the list of all plugin classes currently loaded"]],
        'waitJob': [[PyTango.DevLong, "job id"], [PyTango.DevString, "job state"]],
        }


    #    Attribute definitions
    attr_list = {
        'jobSuccess':
            [[PyTango.DevLong,
            PyTango.SCALAR,
            PyTango.READ]],
        'jobFailure':
            [[PyTango.DevLong,
            PyTango.SCALAR,
            PyTango.READ]],
        'statisticsCollected':
            [[PyTango.DevString,
            PyTango.SCALAR,
            PyTango.READ]],
        'serialize':
            [[PyTango.DevBoolean,
            PyTango.SCALAR,
            PyTango.READ_WRITE]],
    }

    def __init__(self, name):
        PyTango.DeviceClass.__init__(self, name)
        self.set_type(name)
        logger.debug("In DahuDSClass  constructor")
