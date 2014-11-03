#!/usr/bin/env python
# coding: utf-8
from __future__ import with_statement, print_function

__doc__ = """

Data analysis Tango device server ... for UPBL09a

"""
__author__ = "Jérôme Kieffer"
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "23/05/2014"
__status__ = "beta"
__docformat__ = 'restructuredtext'

import sys
import os
import json
import threading
import logging
import time
import types
import multiprocessing
import gc
if sys.version > (3, 0):
    from queue import Queue
else:
    from Queue import Queue
logger = logging.getLogger("dahu.server")
# set loglevel at least at INFO
if logger.getEffectiveLevel() > logging.INFO:
    logger.setLevel(logging.INFO)
logger.setLevel(logging.DEBUG)
import numpy
import PyTango
from .job import Job, plugin_factory

try:
    from rfoo.utils import rconsole
    rconsole.spawn_server()
except ImportError:
    logger.debug("No socket opened for debugging -> please install rfoo")


class DahuDS(PyTango.Device_4Impl):
    """
    Tango device server launcher for Dahu server.
    """
    def __init__(self, cl, name):
        PyTango.Device_4Impl.__init__(self, cl, name)
        self.init_device()
        self.queue = Queue() #queue containing jobs to process
        self.processing_lock = threading.Semaphore()
        self.stat_lock = threading.Semaphore()
        self.last_stats = "No statistics collected yet, please use the 'collectStatistics' method first"
        self.last_failure = -1
        self.last_success = -1
        self.statistics_threads = None
#         self._ncpu_sem = threading.Semaphore(multiprocessing.cpu_count())

    def get_name(self):
        """Returns the name of the class"""
        return self.__class__.__name__

    def delete_device(self):

        logger.debug("[Device delete_device method] for device %s" % self.get_name())

    def init_device(self):
        logger.debug("In %s.init_device()" % self.get_name())

        self.set_state(PyTango.DevState.ON)
        self.get_device_properties(self.get_device_class())
        self.set_change_event("jobSuccess", True, False)
        self.set_change_event("jobFailure", True, False)
        self.set_change_event("statisticsCollected", True, False)

    def always_executed_hook(self):
        pass

    def read_attr_hardware(self, data):
        logger.debug("In %s.read_attr_hardware()" % self.get_name())

    def read_jobSuccess(self, attr):
        attr.set_value(self.last_success)

    def read_jobFailure(self, attr):
        attr.set_value(self.last_failure)

    def read_statisticsCollected(self, attr):
        attr.set_value(self.last_stats)

    def getJobState(self, jobId):
        return Job.getStatusFromID(jobId)

    def cleanJob(self, jobId):
        return Job.cleanJobFromID(jobId)

    def listPlugins(self):
        """
        List all plugin currently loaded .... with a brief description
        """
        logger.debug("In %s.listPlugins" % (self.get_name()))
        res = ["List of all plugin currently loaded (use initPlugin to loaded additional plugins):"]
        plugins = list(plugin_factory.registry.keys())
        plugins.sort()
        return os.linesep.join(res+[" %s : %s"%(i,plugin_factory.registry[i].__doc__.split("\n")[0]) for i in plugins])

    def initPlugin(self, name):
        """
        Creates a job with the given plugin
        """
        logger.debug("In %s.initPlugin(%s)" % (self.get_name(), name))
        err = None
        try:
            plugin = plugin_factory(name)
        except Exception as error:
            err = "plugin %s failed to be instanciated: %s" % (name, error)
            logger.error(err)
        if plugin is None or err:
            return "Plugin not found: %s, err" % (name, err)
        else:
            return "Plugin loaded: %s%s%s" % (name,os.linesep,plugin.__doc__)

    def abort(self, jobId):
        """
        Aborts a job

        @param  jobId: ID of the job to stop
        """
        pass

    def quitDahu(self):
        logger.debug("In %s.quitDahu()" % self.get_name())
        logger.info("Quitting DahuDS")
        sys.exit()

    def startJob(self, argin):
        """
        Starts a job

        @param argin: 2-list [<Dahu plugin to execute>, <JSON serialized dict>]
        @return: jobID which is an int (-1 for error)
        """
        logger.debug("In %s.startJob()" % self.get_name())
        name, data_input = argin[:2]
        print(name, data_input)
        if data_input.strip() == "":
            return -1
        job = Job(name, data_input)
        print(job)
        print(job.input_data)
        print(job)
        if job is None:
            return -1
        self.queue.put(job)
        if self.processing_lock._Semaphore__value > 0 :
            t = threading.Thread(target=self.process_queue)
            t.start()
        return job.id

    def process_queue(self):
        """
        Process all jobs in the queue.
        """
        with self.processing_lock:
            while not self.queue.empty():
#                 self._ncpu_sem.acquire()
                job = self.queue.get()
                job.connect_callback(self.finished_processing)
                job.start()

    def finished_processing(self, job):
        """
        callback: when processing is done

        @param job: instance of dahu.job.Job
        """
        logger.debug("In %s.finished_processing id:%s (%s)" % (self.get_name(), job.id, job.status))
#         self._ncpu_sem.release()
        job.clean(wait=False)
        if job.status == job.STATE_SUCCESS:
            self.last_success = job.id
            self.push_change_event("jobSuccess", job.id)
        else:
            sys.stdout.flush()
            sys.stderr.flush()
            self.last_failure = job.id
            self.push_change_event("jobFailure", job.id)
        self.queue.task_done()
        gc.collect()

#TODO one day
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
            self.last_stats += os.linesep + "Statistics collected on %s, the collect took: %.3fs" % (time.asctime(), time.time() - fStartStat)
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
    }

    def __init__(self, name):
        PyTango.DeviceClass.__init__(self, name)
        self.set_type(name);
        logger.debug("In DahuDSClass  constructor")


