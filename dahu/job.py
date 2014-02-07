#!/usr/bin/env python
# -*- coding: utf8 -*-
#
from __future__ import with_statement, print_function

__doc__ = """Contains the Job class which handles jobs. 
            A static part of the class contains statistics of the class
            """
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "GPLv3+"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20140207"
__status__ = "development"

import threading 
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

class Job(object):
    """
    Class Job
    
    * Each instance will be a job
    * Constructor takes an input data and generates the JobId
    * Each instance will gave a "getOutput" method with optional join 
    * there could be a "join" method, waiting for the job to finish
    * Each instance will have a "execute" method  and  returning a JobId 
    * Each instance will have a "setCallBack" method  that stores the name of the external callback 
    * provide status of a job
    
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
    """
    STATE_UNITIALIZED = "uninitialized"
    STATE_RUNNING = "running"
    STATE_SUCCESS = "success"
    STATE_FAILURE = "failure"

    _dictJobs = {}
    _semaphore = threading.Semaphore()
    _fStartTime = time.time()
    _id_class = 0
    _storage_dir = tempfile.tempdir()

    def __init__(self, input_data):
        """
        Constructor of the class Job
        
        @param input_data: Should be a dictionary or a JSON string representing that dictionary
        """
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
        self._sem = threading.Semaphore()
        self._process = None
        self._runtime = None
        self._start = time.time()

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


    def execute(self):
        """
        Launch the processing
        """
        if not self.__bXmlInputSet:
            logger.warning("Not executing job %s as input is empty" % self._jobId)

        if (self.__edPlugin is not None):
            with self._sem:
                self.__edPlugin.connectSUCCESS(self.successPluginExecution)
                self.__edPlugin.connectFAILURE(self.failurePluginExecution)
                self._status = Job.PLUGIN_STATE_RUNNING
                self.__edPlugin.execute()
                return self._jobId
        else:
            logger.warning("Trying to run a plugin that does not exist: %s " % self.__strPluginName)


    def synchronize(self):
        """
        Synchronize the execution of the job with the calling thread.
        """
        with self._sem:
            strStatus = self._status
        if strStatus == Job.PLUGIN_STATE_RUNNING:
            self.__edPlugin.synchronize()
        elif strStatus == Job.PLUGIN_STATE_UNITIALIZED:
            logger.warning("Unable to synchronize %s jobs" % strStatus)
        else:
            self.DEBUG("Unable to synchronize %s jobs" % strStatus)


    @classmethod
    def synchronizeAll(cls):
        """
        Wait for all jobs to finish.
        """
        logger.debug("Job.synchronizeAll class method ")
        listJob = cls._dictJobs.keys()
        for jobid in listJob:
            job = cls._dictJobs[jobid]
            job.synchronize()
        if len(cls._dictJobs) != len(listJob):
            logger.warning("Job.synchronizeAll: New jobs have been launched while synchronizing")


    def successPluginExecution(self, _edObject=None):
        """
        Method called when the execution of the plugin succeeds 
        """
        with self._sem:
            self._status = Job.PLUGIN_STATE_SUCCESS
            logger.info("Plugin %s: success after %.3fs" % (self._jobId, _edObject.getRunTime()))
        try:
            self.__edSlotSUCCESS.call(self._jobId)
        except Exception:
            logger.error("Error in execution of Success call-back for %s" % self._jobId)
            self.writeErrorTrace()
        try:
            self.__edSlotCallBack.call(self._jobId)
        except Exception:
            logger.error("Error in execution of Common call-back (after success) for %s" % self._jobId)
            self.writeErrorTrace()


    def failurePluginExecution(self, _edObject=None):
        """
        Method called when the execution of the plugin failed 
        """
        with self._sem:
            self._status = Job.PLUGIN_STATE_FAILURE
            logger.info("Plugin %s: failure after %.3fs" % (self._jobId, _edObject.getRunTime()))
        try:
            self.__edSlotFAILURE.call(self._jobId)
        except Exception:
            logger.error("Error in execution of Failure call-back for %s" % self._jobId)
            self.writeErrorTrace()
        try:
            self.__edSlotCallBack.call(self._jobId)
        except Exception:
            logger.error("Error in execution of Common call-back (after failure) for %s" % self._jobId)
            self.writeErrorTrace()


    def connectSUCCESS(self, _oMethod):
        """
        @param _oMethod: function or method to be called - back
        """

        with self._sem:
            if (_oMethod != None):
                self.__edSlotSUCCESS.connect(_oMethod)


    def connectFAILURE(self, _oMethod):
        """
        @param _oMethod: function or method to be called - back
        """
        with self._sem:
            if (_oMethod != None):
                self.__edSlotFAILURE.connect(_oMethod)


    def connectCallBack(self, _oMethod):
        """
        @param _oMethod: function or method to be called - back
        """
        with self._sem:
            if (_oMethod != None):
                self.__edSlotCallBack.connect(_oMethod)

    @property
    def jobId(self):
        """
        @return: JobId 
        @rtype: integer
        """
        return self._jobId

    @property
    def process(self):
        """
        @return: the processing instance
        @rtype: python object
        """
        return self._process

    @property
    def status(self):
        """
        @return: status of the Job
        @rtype: string
        """
        return self._status

    def getName(self):
        return self.__name
    def setName(self, _strName):
        if self.__name is None:
            self.__name = _strName
        else:
            logger.warning("Job.setName: One cannot rename a Job !!!")
    name = property(getName, setName, "Job.name: nickname of the job")




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


    def cleanJob(self, forceGC=True):
        """
        Frees the memory associated with the top level plugin
        @param forceGC: Force garbage collection after clean-up
        @type forceGC: boolean
        """
        self.synchronize()
        with self._sem:
            if self._process is not None:
                self.__pathXSDOutput = self.__edPlugin.strPathDataOutput
                self.__pathXSDInput = self.__edPlugin.strPathDataInput
                self.__runtime = self.__edPlugin.getRunTime()
                self.__edPlugin = None
        if forceGC:
            gc.collect()


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
        lstStrOut = []
        output = []
        fExecTime = time.time() - cls._fStartTime
        keys = cls._dictJobs.keys()
        keys.sort()
        for num, key in enumerate(keys) :
            job = cls._dictJobs[key]
            if job.getPlugin() is None:
                runtime = job.__runtime
            else:
                runtime = job.getPlugin().getRunTime()
            output.append([num, key, job.getStatus(), runtime, job.getMemSize()])
        output.sort()
        iNbJob = max(1, len(keys))
        lstStrOut.append("_" * 110)
        lstStrOut.append("%s\t|\t%s\t\t\t\t|\t%s\t|\t%s\t\t|\t%s" % ("nr", "EDPluginName-Id", "status", "runtime", "memory"))
        lstStrOut.append("_" * 110)
        fWall = 0.0
        fSumProd = 0.0
        fSumX = 0.0
        fSumX2 = 0.0
        for oneJob in output:
            fWall += oneJob[3]
            fSumX += oneJob[0]
            fSumX2 += oneJob[0] * oneJob[0]
            fSumProd += oneJob[0] * oneJob[3]
            lstStrOut.append("%s\t|\t%s\t|\t%s\t|\t%9.3f\t|\t%s" % tuple(oneJob))
        lstStrOut.append("_" * 110)
        lstStrOut.append("Total execution time (Wall): %.3fs, Execution time: %.3fs. SpeedUp: %.3f" % (fWall, fExecTime, fWall / fExecTime))
        lstStrOut.append("Average execution time (Wall/N): %.3fs, Average throughput: %.3fs" % (fWall / iNbJob, fExecTime / iNbJob))
        if len(keys) > 1:
            fSlope = (iNbJob * fSumProd - fSumX * fWall) / (iNbJob * fSumX2 - fSumX * fSumX)
            fOrd = (fWall - fSlope * fSumX) / iNbJob
        else:
            fSlope = 0.0
            fOrd = fWall
        lstStrOut.append("Regression of execution time: ExecTime = %.3f + %f * NbJob" % (fOrd, fSlope))
        strOutput = os.linesep.join(lstStrOut)
        EDVerbose.screen(strOutput)
        return strOutput

