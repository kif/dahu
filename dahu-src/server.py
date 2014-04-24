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
__date__ = "03/03/2014"
__status__ = "beta"
__docformat__ = 'restructuredtext'

import os
import json
import threading
import logging
import time
import types
logger = logging.getLogger("dahu.server")
# set loglevel at least at INFO
if logger.getEffectiveLevel() > logging.INFO:
    logger.setLevel(logging.INFO)
import numpy
import PyTango

class DahuDS(PyTango.Device_4Impl):
    """
    Tango device server launcher for Dahu server.
    """
    def __init__(self, cl, name):
        PyTango.Device_4Impl.__init__(self, cl, name)
        self.init_device()
        self.jobQueue = Queue()
        self.processingSem = threading.Semaphore()
        self.statLock = threading.Semaphore()
        self.lastStatistics = "No statistics collected yet, please use the 'collectStatistics' method first"
        self.lastFailure = "No job Failed (yet)"
        self.lastSuccess = "No job succeeded (yet)"
        self.statistics_threads = threading.Thread()
        self.statistics_threads.join()

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
        attr.set_value(self.lastSuccess)

    def read_jobFailure(self, attr):
        attr.set_value(self.lastFailure)

    def read_statisticsCollected(self, attr):
        attr.set_value(self.lastStatistics)

    def getJobState(self, jobId):
        return Job.getStatusFromID(jobId)

    def cleanJob(self, jobId):
        return Job.cleanJobFromID(jobId)

    def initPlugin(self, strPluginName):
        plugin = EDFactoryPluginStatic.loadPlugin(strPluginName)
        if plugin is None:
            return "Plugin not found: %s" % strPluginName
        else:
            return "Plugin loaded: %s" % strPluginName

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

    def start(self, argin):
        """
        Starts a job
        @param argin: 2-list [ "EDPluginName", "<xml/><XSDataInputPluginName>...."]
        @return: jobID which is a sting: Plugin-000001
        """
        logger.debug("In %s.startJob()" % self.get_name())
        name, xsd = argin[:2]
        if xsd.strip() == "":
            return
        edJob = Job(name)
        if edJob is None:
            return "Error in load Plugin"
        jobId = edJob.getJobId()
        edJob.setDataInput(xsd)
        self.jobQueue.put(edJob)
        if self.processingSem._Semaphore__value > 0 :
            t = threading.Thread(target=self.startProcessing)
            t.start()
        return jobId

    def startProcessing(self):
        """
        Process all jobs in the queue.
        """
        with self.processingSem:
            while not self.jobQueue.empty():
                self.__semaphoreNbThreads.acquire()
                edJob = self.jobQueue.get()
                edJob.connectSUCCESS(self.successJobExecution)
                edJob.connectFAILURE(self.failureJobExecution)
                edJob.execute()

    def successJobExecution(self, jobId):
        logger.debug("In %s.successJobExecution(%s)" % (self.get_name(), jobId))
        with self.locked():
            self.__semaphoreNbThreads.release()
            Job.cleanJobfromID(jobId, False)
            self.lastSuccess = jobId
            self.push_change_event("jobSuccess", jobId)
            gc.collect()

    def failureJobExecution(self, jobId):
        logger.debug("In %s.failureJobExecution(%s)" % (self.get_name(), jobId))
        with self.locked():
            self.__semaphoreNbThreads.release()
            Job.cleanJobfromID(jobId, False)
            self.lastFailure = jobId
            self.push_change_event("jobFailure", jobId)
            sys.stdout.flush()
            sys.stderr.flush()
            gc.collect()

    def getRunning(self):
        """
        retrieve the list of plugins currently under execution (with their plugin-Id)
        """
        return EDStatus.getRunning()

    def getSuccess(self):
        """
        retrieve the list of plugins finished with success (with their plugin-Id)
        """
        return EDStatus.getSuccess()

    def getFailure(self):
        """
        retrieve the list of plugins finished with failure (with their plugin-Id)
        """
        return EDStatus.getFailure()

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
        with self.statLock:
            fStartStat = time.time()
            self.lastStatistics = Job.stats()
            self.lastStatistics += os.linesep + "Statistics collected on %s, the collect took: %.3fs" % (time.asctime(), time.time() - fStartStat)
            self.push_change_event("statisticsCollected", self.lastStatistics)


    def getStatistics(self):
        """
        just return statistics previously calculated
        """
        self.statistics_threads.join()
        return self.lastStatistics

    def getJobOutput(self, jobId):
        """
        Retrieve XML output form a job
        @param jobId: name of the job
        @return: output from a job
        """
        return Job.getDataOutputFromId(jobId)

    def getJobInput(self, jobId):
        """
        Retrieve XML input from a job
        @param jobId: name of the job
        @return: xml input from a job
        """
        return Job.getDataInputFromId(jobId)

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
        'startJob': [[PyTango.DevVarStringArray, "[<Dahu plugin to execute>,<XML input>]"], [PyTango.DevString, "job id"]],
        'abort': [[PyTango.DevString, "job id"], [PyTango.DevBoolean, ""]],
        'getJobState': [[PyTango.DevString, "job id"], [PyTango.DevString, "job state"]],
        "initPlugin": [[PyTango.DevString, "plugin name"], [PyTango.DevString, "Message"]],
        "cleanJob":[[PyTango.DevString, "job id"], [PyTango.DevString, "Message"]],
        "collectStatistics":[[PyTango.DevVoid, "nothing needed"], [PyTango.DevVoid, "Collect some statistics about jobs within Dahu"]],
        "getStatistics":[[PyTango.DevVoid, "nothing needed"], [PyTango.DevString, "Retrieve statistics about Dahu-jobs"]],
        'getJobOutput': [[PyTango.DevString, "job id"], [PyTango.DevString, "job output XML"]],
        'getJobInput': [[PyTango.DevString, "job id"], [PyTango.DevString, "job input XML"]],
        }


    #    Attribute definitions
    attr_list = {
        'jobSuccess':
            [[PyTango.DevString,
            PyTango.SCALAR,
            PyTango.READ]],
        'jobFailure':
            [[PyTango.DevString,
            PyTango.SCALAR,
            PyTango.READ]],
        "statisticsCollected":
            [[PyTango.DevString,
            PyTango.SCALAR,
            PyTango.READ]],
    }

    def __init__(self, name):
        PyTango.DeviceClass.__init__(self, name)
        self.set_type(name);
        logger.debug("In DahuDSClass  constructor")


