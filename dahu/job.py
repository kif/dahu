#!/usr/bin/env python
# -*- coding: utf8 -*-
#
from __future__ import with_statement, print_function

__doc__ = """Contains the Job class which handles jobs. 
            A static part of the class contains statistics of the class
            """
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20140207"
__status__ = "development"

from threading import Thread, Semaphore 
import time
import os
import sys
import gc
import types
import json
import numpy
import tempfile
import logging
logger = logging.getLogger("job")

class Job(Thread):
    """
    Class Job
    
    * Each instance will be a job
    * Constructor takes an input data and generates the JobId
    * Each instance will gave a "getOutput" method with optional join 
    * there could be a "join" method, waiting for the job to finish
    * Each instance will have a "execute" method  and  returning a JobId 
    * Each instance will have a "setCallBack" method  that stores the name of the external callback 
    * provide status of a job
    * Each instance has an abort method which can be used to stop processing (or a server)
    
    
    Static part:
    * keeps track of all jobs status
    * leave the time to job to initialize
    * static class retrieve job-instance, status, small-log ...
    * does not manage workload of the computer, should be managed at the ExecPlugin level
    
    Used for the tango binding
    
    == class variables ==
    dictPluginStatus[pluginName] = ["uninitialized"|"running"|"executed"|"failed"]
    dictJobs [JobId] = Job.Instance
    
    == static methods ==
    getJob(JobId)
    
    
    RESERVED keywords from Thread:
    start, run, join, name, ident, is_alive, daemon
    
    start is overridden with a call to the factory to instanciate the plugin
    
    """
    STATE_UNITIALIZED = "uninitialized"
    STATE_RUNNING = "running"
    STATE_SUCCESS = "success"
    STATE_FAILURE = "failure"
    STATE_ABORTED = "aborted"
    STATE = [STATE_UNITIALIZED, STATE_RUNNING, STATE_SUCCESS, STATE_FAILURE, STATE_ABORTED]

    _dictJobs = {}
    _semaphore = Semaphore()
    _global_start_time = time.time()
    _id_class = 0
    _storage_dir = tempfile.tempdir()

    def __init__(self, input_data):
        """
        Constructor of the class Job
        
        @param input_data: Should be a dictionary or a JSON string representing that dictionary
        """
        Thread.__init__(self)
        if type(input_data) in types.StringTypes:
            if os.path.isfile(input_data):
                self._input_data = json.load(open(input_data))
            else:
                self._input_data = json.loads(input_data)
        else:
            self._input_data = dict(input_data)
        self._status = Job.PLUGIN_STATE_UNITIALIZED
        with self.__class__._semaphore:
            self.__class__._id_class += 1
            self._jobId = self.__class__._id_class
            self.__class__._dictJobs[self.__class__._id_class] = self
        self._input_data["id"] = self._jobId
        self.data_on_disk = False
        self._output_data = None
        self._sem = Semaphore()
        self._plugin = None
        self._runtime = None
        self._start = time.time()
        self._runtime = None
        self._name = self._input_data.get("name", "Plugin")
        # list of methods to be called at the end of the processing
        self._callbacks = []


    @property
    def input_data(self):
        """
        Returns the job input data
        """
        if self.data_on_disk:
            return json.load(open(self.data_on_disk + ".in"))
        else:
            return self._input_data

    @property
    def output_data(self):
        """
        Returns the job output data
        @param _bWait: shall we wait for the plugin to finish to retrieve output data: Yes by default.
        @type _bWait: boolean
        """
        if self._status in [self.STATE_SUCCESS, self.STATE_SUCCESS]:
            if self.data_on_disk:
                return json.load(open(self.data_on_disk + ".out"))
            else:
                return self._output_data
        else:
            logger.warning("Getting output_data for job id %d in state %s." % (self._jobId, self._status))
            return self._output_data

    def start(self):
        """
        We need to create the plugin before starting the new tread... 
        """
        try:
            self._plugin = plugin_factory(self._name)
        except Exception as error:
            self._log_error("plugin %s failed to be instanciated."%self._name)
            self._run_callbacks()
        else:
            #finally launch the new thread.
            Thread.start(self)

    def abort(self):
        """
        Tell the job to stop !
        To be implemented
        """
        with self._sem:
            self._status = self.STATE_ABORTED
            self._run_teardown()
            self._run_callbacks()
        
    def run(self):
        """
        Defines the sequence of execution of the plugin
        1) the the state to "running"
        2) sets the input data to the plugin
        3) run the set-up
        4) run the process
        4) run the tear-down
        5) run the call-backs
        """
        self._status = self.STATE_RUNNING
        self._run_setup()
        if self._status == self.STATE_FAILURE:
            self._run_callbacks()
            return
        self._run_process()
        if self._status == self.STATE_FAILURE:
            self._run_callbacks()
            return
        self._run_teardown()
        if self._status == self.STATE_FAILURE:
            self._run_callbacks()
            return
        self._run_callbacks()
        if self._status == self.STATE_RUNNING:
            self._status = self.STATE_SUCCESS

    def _run_setup(self):
        "TODO"
        setup_name = self._plugin.DEFAULT_SET_UP

        for cb in self._callbacks:
            if "__call__" in dir(cb):
                try:
                    cb(self)
                except Exception as error:
                    self._log_error("Error while calling %s" % cb)
    def _run_process(self):
        "TODO"
        process_name = self._plugin.DEFAULT_TEAR_DOWN

        for cb in self._callbacks:
            if "__call__" in dir(cb):
                try:
                    cb(self)
                except Exception as error:
                    self._log_error("Error while calling %s" % cb)
    def _run_teardown(self):
        "TODO"
        teardown_name = self._plugin.DEFAULT_TEAR_DOWN

        for cb in self._callbacks:
            if "__call__" in dir(cb):
                try:
                    cb(self)
                except Exception as error:
                    self._log_error("Error while calling %s" % cb)
    def _run_callbacks(self):
        for cb in self._callbacks:
            if "__call__" in dir(cb):
                try:
                    cb(self)
                except Exception as error:
                    self._log_error("Error while calling %s"%cb)
        self._status = self.STATE_SUCCESS
                    
    def _log_error(self, msg):
        """
        log an error message in the output 
        """
        exc_type, exc_value, exc_traceback = sys.exc_info()
        err_msg = [msg, "%s: %s" % (exc_type, exc_value)]
        for line in traceback.extract_tb(exc_traceback):
            err_msg.append("  File \"%s\", line %d, in %s" % (line[0], line[1], line[2]))
            err_msg.append("\t\t%s" % line[3])
        with self._sem:
            self._status = self.STATE_FAILURE
            if self._output_data is None:
                self._output_data = {}
            if "error" not in self._output_data:
                self._output_data["error"] = err_msg
            else:
                self._output_data["error"] += ["*"*50] + err_msg

#
#    def execute(self):
#        """
#        Launch the processing
#        """
#        if not self.__bXmlInputSet:
#            logger.warning("Not executing job %s as input is empty" % self._jobId)
#
#        if (self.__edPlugin is not None):
#            with self._sem:
#                self.__edPlugin.connectSUCCESS(self.successPluginExecution)
#                self.__edPlugin.connectFAILURE(self.failurePluginExecution)
#                self._status = Job.PLUGIN_STATE_RUNNING
#                self.__edPlugin.execute()
#                return self._jobId
#        else:
#            logger.warning("Trying to run a plugin that does not exist: %s " % self.__strPluginName)


#    def synchronize(self, timeout=None):
#        """
#        Synchronize the execution of the job with the calling thread.
#        """
#        with self._sem:
#            strStatus = self._status
#        if strStatus == Job.PLUGIN_STATE_RUNNING:
#            self.__edPlugin.synchronize()
#        elif strStatus == Job.PLUGIN_STATE_UNITIALIZED:
#            logger.warning("Unable to synchronize %s jobs" % strStatus)
#        else:
#            self.DEBUG("Unable to synchronize %s jobs" % strStatus)



    def connect_callback(self, method=None):
        """
        @param method: function or method to be called - back
        """
        if method:
            with self._sem:
                if "__call__" in dir(method):
                    self._callbacks.append(method)
                else:
                    logger.error("Non callable callback method: %s" % method)


    def cleanJob(self, force=True):
        """
        Frees the memory associated with the plugin
        
        @param force: Force garbage collection after clean-up
        TODO
        """
        self.synchronize()
        with self._sem:
            if self._plugin is not None:
                self.__pathXSDOutput = self.__edPlugin.strPathDataOutput
                self.__pathXSDInput = self.__edPlugin.strPathDataInput
                self.__runtime = self.__edPlugin.getRunTime()
                self.__edPlugin = None
        if force:
            gc.collect()



################################################################################
# Properties
################################################################################
    @property
    def jobId(self):
        """
        @return: JobId 
        @rtype: integer
        """
        return self._jobId

    @property
    def plugin(self):
        """
        @return: the processing instance
        @rtype: python object
        """
        return self._plugin

    @property
    def status(self):
        """
        @return: status of the Job
        @rtype: string
        """
        return self._status

    def getName(self):
        return self._name
    def setName(self, name):
        if self._name is None:
            self._name = name
        else:
            logger.error("Job.setName: One cannot rename a Job !!!")
    name = property(getName, setName, "Job.name: nickname of the job")


################################################################################
# Class methods
################################################################################
    @classmethod
    def synchronize_all(cls):
        """
        Wait for all jobs to finish.
        """
        logger.debug("Job.synchronize_all class method ")
        listJob = cls._dictJobs.keys()
        for jobid in listJob:
            job = cls._dictJobs[jobid]
            job.synchronize()
        if len(cls._dictJobs) != len(listJob):
            logger.warning("Job.synchronize_all: New jobs have been launched while synchronizing")

    @classmethod
    def getStatusFromID(cls, jobId):
        """
        Retrieve the job (hence the plugin) status
        
        @param jobId: the Job identification number
        @type jobId: string
        @return: the Job status 
        @rtype: string 
        """
        if jobId in cls._dictJobs:
            strRet = cls._dictJobs[jobId].getStatus()
        else:
            strRet = "Unable to retrieve such job: %s" % jobId
            logger.warning(strRet)
        return strRet
    getStatusFromId = getStatusFromID


    @classmethod
    def getJobFromID(cls, jobId):
        """
        Retrieve the job (hence the plugin)
        
        @param jobId: the Job identification number
        @return: the "Job instance", which contains the plugin and the status
        @rtype: a Python object, instance of Job. 
        """
        if jobId in cls._dictJobs:
            return cls._dictJobs[jobId]
        else:
            logger.warning("Unable to retrieve such Job: %s" % jobId)
    getJobFromId = getJobFromID



    @classmethod
    def cleanJobfromId(cls, jobId, forceGC=True):
        """
        Frees the memory associated with the top level plugin
        
        @param jobId: the Job identification number
        @type jobId: string
        @param forceGC: Force garbage collection after clean-up
        @type forceGC: boolean
        """
        if jobId in cls._dictJobs:
            job = cls._dictJobs[jobId]
            job.cleanJob(forceGC)
            strRet = "Job %s cleaned" % jobId
        else:
            strRet = "Unable to retrieve such Job: %s" % jobId
            logger.warning(strRet)
        return strRet
    cleanJobfromID = cleanJobfromId


    @classmethod
    def getDataOutputFromId(cls, jobId):
        """
        Returns the Plugin Output Data
        @param jobId: job idenfier 
        @type jobId: string
        @return: Job.DataOutput XML string
        """
        output = None
        job = cls.getJobFromId(jobId)
        if job is not None:
            output = job.getDataOutput()
        return output or ""
    getDataOutputFromID = getDataOutputFromId


    @classmethod
    def getDataInputFromId(cls, jobId):
        """
        Returns the Plugin Input Data
        @param jobId: job idenfier 
        @type jobId: string
        @return: Job.DataInput XML string
        """
        output = None
        job = cls.getJobFromId(jobId)
        if job is not None:
            output = job.getDataInput()
        return output or ""
    getDataInputFromID = getDataInputFromId


    @classmethod
    def stats(cls):
        """
        Retrieve some statistics and print them
        """
        lstStrOut = [""]
        output = []
        run_time = time.time() - cls._global_start_time
        keys = cls._dictJobs.keys()
        keys.sort()

        output = [ (i, cls._dictJobs[i]._name, cls._dictJobs[i]._status, cls._dictJobs[i]._runtime)
                  for i, k in enumerate(keys)]
        total_jobs = max(1, len(keys))
        lstStrOut.append("_" * 110)
        lstStrOut.append("%s\t|\t%s\t\t\t\t|\t%s\t|\t%s" % ("id", "Name", "Status", "run time"))
        lstStrOut.append("_" * 110)
        wall_time = 0.0
        fSumProd = 0.0
        fSumX = 0.0
        fSumX2 = 0.0
        for oneJob in output:
            wall_time += oneJob[3]
            fSumX += oneJob[0]
            fSumX2 += oneJob[0] * oneJob[0]
            fSumProd += oneJob[0] * oneJob[3]
            lstStrOut.append("%s\t|\t%s\t|\t%s\t|\t%9.3f\t|\t%s" % tuple(oneJob))
        lstStrOut.append("_" * 110)
        lstStrOut.append("Total execution time (Wall): %.3fs, Execution time: %.3fs. SpeedUp: %.3f" % (wall_time, run_time, wall_time / run_time))
        lstStrOut.append("Average execution time (Wall/N): %.3fs, Average throughput: %.3fs" % (wall_time / total_jobs, run_time / total_jobs))
        if len(keys) > 1:
            fSlope = (total_jobs * fSumProd - fSumX * wall_time) / (iNbJob * fSumX2 - fSumX * fSumX)
            fOrd = (wall_time - fSlope * fSumX) / iNbJob
        else:
            fSlope = 0.0
            fOrd = wall_time
        lstStrOut.append("Regression of execution time: ExecTime = %.3f + %f * NbJob" % (fOrd, fSlope))
        strOutput = os.linesep.join(lstStrOut)
        logger.info(strOutput)
        return strOutput

