#!/usr/bin/env python3
#
"""
Data Analysis RPC server over Tango:

Contains the Job class which handles jobs.
A static part of the class contains statistics of the class
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "11/03/2026"
__status__ = "production"

import gc
import json
import logging
import os
import sys
import time
import traceback
from threading import Semaphore, Thread

import six

from . import utils
from .factory import plugin_factory
from .utils import NumpyEncoder

logger = logging.getLogger("dahu.job")
# logger.setLevel(logging.DEBUG)

# Python 2to3 compatibility
StringTypes = (six.binary_type, six.text_type)


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
    STATE_INEXISTANT = "inexistant"
    STATE_UNINITIALIZED = "uninitialized"
    STATE_STARTING = "starting"
    STATE_RUNNING = "running"
    STATE_SUCCESS = "success"
    STATE_FAILURE = "failure"
    STATE_ABORTED = "aborted"
    STATE = [STATE_UNINITIALIZED, STATE_STARTING, STATE_RUNNING, STATE_SUCCESS, STATE_FAILURE, STATE_ABORTED]

    _dictJobs = {}
    _semaphore = Semaphore()
    _global_start_time = time.time()
    _id_class = 0
    _storage_dir = utils.get_workdir()
    if not os.path.isdir(_storage_dir):
        os.makedirs(_storage_dir)

    def __init__(self, name="plugin.Plugin", input_data=None):
        """
        Constructor of the class Job

        :param name: name of the plugin to be instanciated
        :param input_data: Should be a dictionary or a JSON string representing that dictionary
        """
        if input_data is None:
            input_data = {}
        Thread.__init__(self)
        self._status = Job.STATE_UNINITIALIZED
        self._name = name
        if isinstance(input_data, StringTypes):
            if os.path.isfile(input_data):
                with open(input_data) as f:
                    self._input_data = json.load(f)
            else:
                self._input_data = json.loads(input_data)
        else:
            self._input_data = dict(input_data)
        self._input_data["plugin_name"] = self._name
        with self.__class__._semaphore:
            self.__class__._id_class += 1
            self._jobId = self.__class__._id_class
            self.__class__._dictJobs[self.__class__._id_class] = self
        self._input_data["job_id"] = self._jobId
        self.data_on_disk = False
        self._output_data = {}
        self._sem = Semaphore()
        self._plugin = None
        self._runtime = None
        self._start_time = time.time()
        # list of methods to be called at the end of the processing
        self._callbacks = []

    def __repr__(self):
        if self._plugin is None:
            txt = "dahu job (thread) #%i using finished plugin (%s) %s" % (
                  self._jobId, self._name, self._status)

        else:
            txt = "dahu job (thread) #%i using plugin %s currently %s" % (
                  self._jobId, self._plugin.__class__.__name__, self._status)
        return txt

    def start(self):
        """
        We need to create the plugin before starting the new tread...
        """
        self._status = self.STATE_STARTING
        try:
            self._plugin = plugin_factory(self._name)
        except Exception as error:
            self._log_error(f"plugin {self._name} failed to be instanciated, raised: {error}")
            self._run_callbacks()
        else:
            if self._plugin is None:
                self._log_error(f"plugin {self._name} failed to be instanciated.")
                logger.debug(plugin_factory.registry)
                self._run_callbacks()
            else:
                # finally launch the new thread.
                Thread.start(self)

    def join(self, timeout=None):
        if self._status in (self.STATE_RUNNING, self.STATE_STARTING):
            Thread.join(self, timeout)

    def abort(self):
        """
        Tell the job to stop !

        Needs to be implemented into the plugin !

        :return: True if the job was running and has been asked to stop
        """
        if self._status != self.STATE_RUNNING:
            logger.warning(f"Job {self._jobId} is {self._status}, not aborting it")
            return False
        with self._sem:
            self._status = self.STATE_ABORTED
            self._output_data[self._status] = utils.get_isotime()
        # outside of the semaphore: the plugin may fail and log the error
        self._run_("abort")
        return True

    def run(self):
        """
        Defines the sequence of execution of the plugin
        1) the the state to "running"
        2) sets the input data to the plugin
        3) run the set-up
        4) run the process
        4) run the tear-down: always runs tear-down !
        5) run the call-backs
        """
        self._status = self.STATE_RUNNING
        self._plugin.input.update(self._input_data)
        self._run_("setup")
        if self._status != self.STATE_FAILURE:
            self._run_("process")
        self._run_("teardown")
        if self._status != self.STATE_FAILURE:
            self._output_data.update(self._plugin.output)
            if self._status == self.STATE_RUNNING:
                self._status = self.STATE_SUCCESS
        self._run_callbacks()

    def _run_(self, what):
        """
        run setup, process, teardown or abort ...

        :param what: setup, process or teardown
        @parma args: argument list to be passed to the method

        """
        methods = {"process": self._plugin.DEFAULT_PROCESS,
                   "setup": self._plugin.DEFAULT_SET_UP,
                   "teardown": self._plugin.DEFAULT_TEAR_DOWN,
                   "abort": self._plugin.DEFAULT_ABORT    }
        assert what in methods
        name = methods.get(what)
        if name in dir(self._plugin):
            method = self._plugin.__getattribute__(name)
            if "__call__" in dir(method):
#                 if what in self._input_data:
#                     try:
#                         method(self._input_data[what])
#                     except Exception as error:
#                         self._log_error("Error %s while calling %s.%s with argument %s" %
#                                         (error, self._plugin.__class__.__name__,
#                                          what, self._input_data[what]))
#                 else:
                try:
                    method()
                except Exception as error:
                    import traceback
                    err_msg = [traceback.format_exc(limit=10), (""
                               f"Error {error} while calling {self._plugin.__class__.__name__}.{what}")]
                    self._log_error(os.linesep.join(err_msg))
        else:
            logger.error(f"No such method {what} in class {self._plugin.__class__.__name__}")

    def _run_callbacks(self):
        self._update_runtime()
        for cb in self._callbacks:
            if "__call__" in dir(cb):
                try:
                    cb(self)
                except Exception as error:
                    self._log_error(f"Error while calling {cb}: {error}")

    def _log_error(self, msg):
        """
        log an error message in the output
        """
        exc_type, exc_value, exc_traceback = sys.exc_info()
        err_msg = [msg, f"{exc_type}: {exc_value}"]
        for line in traceback.extract_tb(exc_traceback):
            err_msg.append("  File \"%s\", line %d, in %s" % (line[0], line[1], line[2]))
            err_msg.append(f"\t\t{line[3]}")
        with self._sem:
            self._status = self.STATE_FAILURE
            if "error" not in self._output_data:
                self._output_data["error"] = err_msg
            else:
                self._output_data["error"] += ["*" * 50] + err_msg
            logger.error(err_msg)

    def _update_runtime(self):
        with self._sem:
            self._runtime = time.time() - self._start_time

    def connect_callback(self, method=None):
        """
        :param method: function or method to be called - back
        """
        if method:
            with self._sem:
                if "__call__" in dir(method):
                    self._callbacks.append(method)
                else:
                    logger.error(f"Non callable callback method: {method}")

    def clean(self, force=False, wait=True):
        """
        Frees the memory associated with the plugin

        :param force: Force garbage collection after clean-up
        :param wait: wait for job to be finished

        """
        logger.debug(f"In clean {self._plugin}")
        if wait and self.is_alive():
            self.join()
        if self._plugin is not None:
            self._update_runtime()
            with self._sem:
                if self._plugin is not None:
                    if self._plugin.input:
                        self._input_data.update(self._plugin.input)
                    if self._plugin.output:
                        self._output_data.update(self._plugin.output)
                    self._output_data["job_runtime"] = self._runtime
                    dirname = os.path.join(utils.get_workdir(), "%04d" % (self._jobId // 1000))
                    with self._semaphore:
                        if not os.path.exists(dirname):
                            os.mkdir(dirname)
                    base_path = os.path.join(dirname, "%05d_%s" % (self._jobId, self._name))
                    with open(base_path + ".inp", "w") as infile:
                        json.dump(self._input_data, infile, indent=4, cls=NumpyEncoder)
                    with open(base_path + ".out", "w") as infile:
                        json.dump(self._output_data, infile, indent=4, cls=NumpyEncoder)
                    self._plugin = None
                    self.data_on_disk = base_path
                    self._output_data = None
                    self._input_data = None
                    logger.info("Finished job #%s, %s with %s in %.3fs", self._jobId, self._name, self._status, self._runtime)
        if force:
            gc.collect()

    synchronize = join

################################################################################
# Properties
################################################################################
    @property
    def id(self):
        """
        :return: JobId
        @rtype: integer
        """
        return self._jobId

    @property
    def plugin(self):
        """
        :return: the processing instance
        @rtype: python object
        """
        return self._plugin

    @property
    def status(self):
        """
        :return: status of the Job
        @rtype: string
        """
        return self._status

    @property
    def input_data(self):
        """
        Returns the job input data
        """
        with self._sem:
            if self.data_on_disk:
                with open(self.data_on_disk + ".inp") as fp:
                    return json.load(fp)
            else:
                return self._input_data

    @property
    def output_data(self):
        """
        Returns the job output data
        :param _bWait: shall we wait for the plugin to finish to retrieve output data: Yes by default.
        :type _bWait: boolean
        """
        with self._sem:
            if self._status in [self.STATE_SUCCESS, self.STATE_FAILURE, self.STATE_ABORTED]:
                if self.data_on_disk:
                    with open(self.data_on_disk + ".out") as fp:
                        return json.load(fp)
                else:
                    return self._output_data
            else:
                logger.warning("Getting output_data for job id %d in state %s." % (self._jobId, self._status))
                return self._output_data

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
        listJob = list(cls._dictJobs.keys())
        for jobid in listJob:
            job = cls._dictJobs[jobid]
            job.join()
        if len(cls._dictJobs) != len(listJob):
            logger.warning("Job.synchronize_all: New jobs have been launched while synchronizing")

    @classmethod
    def synchronize_job(cls, jobId, timeout=None):
        """
        Wait for all a specific jobs to finish.

        :param jobId: identifier of the job ... intg
        :param timeout: timeout in second to wait
        :return: status of the job
        """
        logger.debug("Job.synchronize_job class method for id=%s (timeout=%s)", jobId, timeout)
        job = cls.getJobFromID(jobId)

        if job is None:
            res = cls.STATE_INEXISTANT
        else:
            job.join(timeout)
            res = job.status
        logger.debug("Job.synchronize_job(jobid=%s) ==>  %s", jobId, res)
        return res

    @classmethod
    def getStatusFromID(cls, jobId):
        """
        Retrieve the job (hence the plugin) status

        :param jobId: the Job identification number
        :type jobId: int
        :return: the Job status
        @rtype: string
        """
        if jobId < 0:
            jobId = len(cls._dictJobs) + jobId + 1
        if jobId in cls._dictJobs:
            strRet = cls._dictJobs[jobId]._status
        else:
            strRet = f"Unable to retrieve such job: {jobId}"
            logger.warning(strRet)
        return strRet

    getStatusFromId = getStatusFromID

    @classmethod
    def getJobFromID(cls, jobId):
        """
        Retrieve the job (hence the plugin)

        :param jobId: the Job identification number
        :return: the "Job instance", which contains the plugin and the status
        @rtype: a Python object, instance of Job.
        """
        if jobId < 0:
            jobId = len(cls._dictJobs) + jobId + 1
        if jobId in cls._dictJobs:
            return cls._dictJobs[jobId]
        else:
            logger.warning(f"Unable to retrieve such Job: {jobId}")

    getJobFromId = getJobFromID

    @classmethod
    def abort_job_from_id(cls, jobId):
        """
        Ask a running job to stop

        The plugin has to honour it: it is expected to check `is_aborted`.

        :param jobId: the Job identification number
        :type jobId: int
        :return: True if the job was running and has been asked to stop
        """
        job = cls.getJobFromID(jobId)
        if job is None:
            return False
        return job.abort()

    @classmethod
    def clean_job_from_id(cls, jobId, forceGC=True):
        """
        Frees the memory associated with the top level plugin

        :param jobId: the Job identification number
        :type jobId: int
        :param forceGC: Force garbage collection after clean-up
        :type forceGC: boolean
        """
        if jobId < 0:
            jobId = len(cls._dictJobs) + jobId + 1
        if jobId in cls._dictJobs:
            job = cls._dictJobs[jobId]
            job.clean(forceGC)
            strRet = f"Job {jobId} cleaned"
        else:
            strRet = f"Unable to retrieve such Job: {jobId}"
            logger.warning(strRet)
        return strRet

    # Deprecated aliases, kept for backward compatibility with existing code:
    cleanJobfromId = clean_job_from_id
    cleanJobfromID = clean_job_from_id
    cleanJobFromId = clean_job_from_id
    cleanJobFromID = clean_job_from_id

    @classmethod
    def getDataOutputFromId(cls, jobId, as_JSON=False):
        """
        Returns the Plugin Output Data
        :param jobId: job idenfier
        :type jobId: int
        :return: Job.DataOutput JSON string
        """
        none = "" if as_JSON else {}
        output = None
        if jobId < 0:
            jobId = len(cls._dictJobs) + jobId + 1
        if jobId in cls._dictJobs:
            job = cls._dictJobs[jobId]
            if job is not None:
                with job._sem:
                    if job.data_on_disk:
                        with open(job.data_on_disk + ".out") as fp:
                            data = fp.read()
                        if as_JSON:
                            output = data
                        else:
                            output = json.loads(data)
                    else:
                        if as_JSON:
                            output = json.dumps(job._output_data, skipkeys=True, allow_nan=True,
                                                indent=None, separators=(',', ':'),
                                                cls=NumpyEncoder)
                        else:
                            output = job._output_data
        else:
            output = f"No such job: {jobId}"
        return output or none

    getDataOutputFromID = getDataOutputFromId

    @classmethod
    def getDataInputFromId(cls, jobId, as_JSON=False):
        """
        Returns the Plugin Input Data
        :param jobId: job idenfier
        :type jobId: int
        :return: Job.DataInput JSON string
        """
        output = None
        none = "" if as_JSON else {}
        if jobId < 0:
            jobId = len(cls._dictJobs) + jobId + 1
        if jobId in cls._dictJobs:
            job = cls._dictJobs[jobId]
            if job is not None:
                with job._sem:
                    if job.data_on_disk:
                        with open(job.data_on_disk + ".inp") as fp:
                            data = fp.read()
                        if as_JSON:
                            output = data
                        else:
                            output = json.loads(data)
                    else:
                        if as_JSON:
                            output = json.dumps(job._input_data, skipkeys=True, allow_nan=True,
                                                indent=None, separators=(',', ':'),
                                                cls=NumpyEncoder)
                        else:
                            output = job._input_data
        else:
            output = f"No such job: {jobId}"
        return output or none

    getDataInputFromID = getDataInputFromId

    @classmethod
    def getErrorFromId(cls, jobId):
        """
        Returns the error messages from plugin
        :param jobId: job idenfier
        :type jobId: int
        :return: error message as a string
        """
        out = cls.getDataOutputFromId(jobId)
        return os.linesep.join(out.get("error", [])).encode("UTF-8")

    getErrorFromID = getErrorFromId

    @classmethod
    def stats(cls):
        """
        Retrieve some statistics and print them
        """
        lout = [""]
        run_time = time.time() - cls._global_start_time
        keys = list(cls._dictJobs.keys())
        keys.sort()
        output = [(k, cls._dictJobs[k]._name, cls._dictJobs[k]._status, cls._dictJobs[k]._runtime)
                  for k in keys]
        total_jobs = max(1, len(keys))
        lout.append("_" * 80)
        lout.append("{}\t|\t{}\t\t\t|\t{}\t|\t{} (sec)".format("Id", "Name", "Status", "Run-time"))
        lout.append("_" * 80)
        wall_time = 0.0
        sum_xy = 0.0
        sum_x = 0.0
        sum_xx = 0.0
        for ajob in output:
            if ajob[3]:
                wall_time += ajob[3]
                sum_x += ajob[0]
                sum_xx += ajob[0] * ajob[0]
                sum_xy += ajob[0] * ajob[3]
            lout.append("{}\t|\t{}\t|\t{}\t|\t{}".format(*tuple(ajob)))
        lout.append("_" * 80)
        lout.append(f"Total execution time (Wall): {wall_time:.3f}s, Execution time: {run_time:.3f}s. SpeedUp: {wall_time / run_time:.3f}")
        lout.append(f"Average execution time (Wall/N): {wall_time / total_jobs:.3f}s, Average throughput: {run_time / total_jobs:.3f}s")
        if len(keys) > 1:
            slope = (total_jobs * sum_xy - sum_x * wall_time) / (len(keys) * sum_xx - sum_x * sum_x)
            ord0 = (wall_time - slope * sum_x) / len(keys)
        else:
            slope = 0.0
            ord0 = wall_time
        lout.append(f"Regression of execution time: ExecTime = {ord0:.3f} + {slope:f} * NbJob")
        sout = os.linesep.join(lout)
        logger.info(sout)
        return sout
