#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""
Data Analysis RPC server over Tango: 

Definiton of plugins
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "01/09/2020"
__status__ = "production"

import os
import logging
import cProfile
import time
from .factory import plugin_factory
from .utils import get_workdir

logger = logging.getLogger("dahu.plugin")


class Plugin(object):
    """
    A plugin is instanciated

    * Gets its input parameters as a dictionary from the setup method
    * Performs some work in the process
    * Sets the result as output attribute, should be a dictionary
    * The process can be an infinite loop or a server which can be aborted using the abort method

    """
    DEFAULT_SET_UP = "setup"  # name of the method used to set-up the plugin (close connection, files)
    DEFAULT_PROCESS = "process"  # specify how to run the default processing
    DEFAULT_TEAR_DOWN = "teardown"  # name of the method used to tear-down the plugin (close connection, files)
    DEFAULT_ABORT = "abort"  # name of the method used to abort the plugin (if any. Tear_Down will be called)
    REPROCESS_IGNORE = [] #list of keys from input to be ignored when reprocessing data 


    def __init__(self):
        """
        We assume an empty constructor
        """
        self.input = {}
        self.output = {}
        self._logging = []  # stores the logging information to send back
        self.is_aborted = False
        self.__profiler = None

    def get_name(self):
        return self.__class__.__name__

    def setup(self, kwargs=None):
        """
        This is the second constructor to setup
        input variables and possibly initialize
        some objects
        """
        if kwargs is not None:
            self.input.update(kwargs)
        if self.input.get("do_profiling"):
            self.__profiler = cProfile.Profile()
            self.__profiler.enable()

    def process(self):
        """
        main processing of the plugin
        """
        pass

    def teardown(self):
        """
        method used to tear-down the plugin (close connection, files)

        This is always run, even if process fails
        """
        self.output["logging"] = self._logging
        if self.input.get("do_profiling"):
            self.__profiler.disable()
            name = "%05i_%s.%s.profile" % (self.input.get("job_id", 0), self.__class__.__module__, self.__class__.__name__)
            profile_file = os.path.join(get_workdir(), name)
            self.log_error("Profiling information in %s" % profile_file, do_raise=False)
            self.__profiler.dump_stats(profile_file)

    def get_info(self):
        """
        """
        return os.linesep.join(self._logging)

    def abort(self):
        """
        Method called to stop a server process
        """
        self.is_aborted = True

    def log_error(self, txt, do_raise=True):
        """
        Way to log errors and raise error
        """
        if do_raise:
            err = "ERROR in %s: %s" % (self.get_name(), txt)
            logger.error(err)
        else:
            err = "Warning in %s: %s" % (self.get_name(), txt)
            logger.warning(err)
        self._logging.append(err)
        if do_raise:
            raise RuntimeError(err)

    def log_warning(self, txt):
        """
        Way to log warning
        """
        err = "Warning in %s: %s" % (self.get_name(), txt)
        logger.warning(err)
        self._logging.append(err)

    def wait_for(self, job_id):
        """Wait for another job to be finished ...

        :param job_id: identifier for the job
        :return: the job object
        """
        from .job import Job
        assert isinstance(job_id, int)
        TIMEOUT = getattr(self, 'TIMEOUT', 10.0) #default timeout to 10s
        if job_id>Job._id_class:
            self.log_warning("Not synchronizing job, invalid id: %s" % job_id)
        else:
            status = Job.synchronize_job(job_id, TIMEOUT)
            abort_time = time.time() + TIMEOUT
            while status == Job.STATE_UNINITIALIZED:
                # Wait for job to start
                time.sleep(0.1)
                status = Job.synchronize_job(job_id, TIMEOUT)
                if time.time() > abort_time:
                    self.log_error("Timeout while waiting other job to finish")
                    break
            if status != Job.STATE_SUCCESS:
                self.log_error("Other job ended in %s: aborting myself" % status)
            return Job.getJobFromId(job_id)


class PluginFromFunction(Plugin):
    """
    Template class to build  a plugin from a function
    """

    def __init__(self):
        """
        :param funct: function to be wrapped
        """
        Plugin.__init__(self)

    def __call__(self, **kwargs):
        """
        Behaves like a normal function: for debugging
        """
        self.input.update(kwargs)
        self.process()
        self.teardown()
        return self.output["result"]

    def process(self):
        if self.input is None:
            logger.warning("PluginFromFunction.process: self.input is None !!! %s", self.input)
        else:
            funct_input = self.input.copy()
            if "job_id" in funct_input:
                funct_input.pop("job_id")
            if "plugin_name" in funct_input:
                funct_input.pop("plugin_name")
            self.output["result"] = self.function(**funct_input)


def plugin_from_function(function):
    """
    Create a plugin class from a given function and registers it into the

    :param function: any function
    :return: plugin name to be used by the plugin_factory to get an instance
    """
    logger.debug("creating plugin from function %s" % function.__name__)
    class_name = function.__module__ + "." + function.__name__
    klass = type(class_name, (PluginFromFunction,),
                 {'function': staticmethod(function),
                 "__doc__": function.__doc__})
    plugin_factory.register(klass, class_name)
    return class_name

