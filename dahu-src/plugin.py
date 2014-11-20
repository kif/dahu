#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
from __future__ import with_statement, print_function, division


__doc__ = """
Data Analysis Highly tailored fror Upbl09a
"""
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "19/11/2014"
__status__ = "development"
version = "0.2.0"
from .factory import plugin_factory, register
from .utils import fully_qualified_name, get_workdir
import os
import logging
import cProfile
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
            profile_file = os.path.join(utils.get_workdir(), "%05i_%s" % (self.input.get("job_id"), fully_qualified_name(self)))
            self.__profiler.dump_stats(os.path.join(get_workdir, profile_file))

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
            err = "ERROR in %s: %s" %( self.get_name(), txt)
            logger.error(err)
        else:
            err = "Warning in %s: %s" %( self.get_name(), txt)
            logger.warning(err)
        self._logging.append(err)
        if do_raise:
            raise RuntimeError(err)



class PluginFromFunction(Plugin):
    """
    Template class to build  a plugin from a function
    """
    def __init__(self):
        """
        @param funct: function to be wrapped
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
            print("PluginFromFunction.process: self.input is None !!! %s", self.input)
        else:
            funct_input =  self.input.copy()
            if "job_id" in funct_input:
                funct_input.pop("job_id")
            if "plugin_name" in funct_input:
                funct_input.pop("plugin_name")
            self.output["result"] = self.function(**funct_input)


def plugin_from_function(function):
    """
    Create a plugin class from a given function and registers it into the

    @param function: any function
    @return: plugin name to be used by the plugin_factory to get an instance
    """
    logger.debug("creating plugin from function %s" % function.__name__)
    class_name = function.__module__ + "." + function.__name__
    klass = type(class_name, (PluginFromFunction,),
                 {'function': staticmethod(function),
                 "__doc__": function.__doc__})
    plugin_factory.register(klass, class_name)
    return class_name


