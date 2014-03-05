#!/usr/bin/env python
# -*- coding: utf8 -*-
#

"""
Data Analysis Highly tailored fror Upbl09a 
"""
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20140304"
__status__ = "development"
version = "0.1"
from __future__ import with_statement, print_function
import os

class Plugin(object):
    """
    A plugin is instanciated
    
    * Gets its input parameters as a dictionary from the setup method
    * Performs some work in the process
    * Sets the result as output attribute, should be a dictionary
    * The process can be an infinite loop or a server which can be aborted using the abort method 
    
    """
    IS_DAHU_PLUGIN = True
    DEFAULT_SET_UP = "setup"      # name of the method used to set-up the plugin (close connection, files)
    DEFAULT_PROCESS = "process"   # specify how to run the default processing
    DEFAULT_TEAR_DOWN = "teardown"# name of the method used to tear-down the plugin (close connection, files)
    DEFAULT_ABORT = "abort"       # name of the method used to abort the plugin (if any. Tear_Down will be called)

    def __init__(self):
        """         
        We assume an empty constructor
        """
        self.input = None
        self.output = {}
        self._logging = [] # stores the logging information to send back
        self.is_aborted = False

    def get_name(self):
        return self.__class__.__name__

    def setup(self, kwargs=None):
        """
        This is the second constructor to setup 
        input variables and possibly initialize
        some objects 
        """
        self.input = kwargs

    def process(self, kargs=None):
        """
        main processing of the plugin
        """
        pass

    def teardown(self):
        """
        method used to tear-down the plugin (close connection, files)
        """
        self.output["logging"] = self._logging

    def get_info(self):
        """
        """
        return os.path.linesep.join(self._logging)

    def abort(self):
        """
        Method called to stop a server process
        """
        self.is_aborted = True

class PluginFunction(Plugin):
    """
    Template class to build  a plugin from a function
    """
    def __init__(self, func):
        """
        @param funct: function to be wrapped  
        """
        Plugin.__init__(self)
        self.function = func

    def get_name(self):
        return "Plugin_from_%s" % self.function.__name__

    def __call__(self, **kwargs):
        """
        Behaves like a normal function: for debugging 
        """
        self.setup(kwargs)
        self.process()
#        self.teardown()
        return self.output

    def process(self):
        self.output = self.function(**self.input)

def plugin_from_function(function):
    """
    Instanciate a plugin from a given function
    """
    inst = PluginFunct(function)
    return inst


if __name__ == "__main__":
    #here I should explain how to run the plugin as stand alone:
    p = Plugin()
    p.setup()
    p.process()
    p.teardown()
