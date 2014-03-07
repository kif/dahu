#!/usr/bin/env python
# -*- coding: utf8 -*-
#

"""
Data Analysis Highly tailored for Upbl09a 
"""
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20140304"
__status__ = "development"
version = "0.1"
from __future__ import with_statement, print_function
import os, sys
from threading import Semaphore
from .utils import get_workdir 


class Factory(object):
    """
    This is a factory, it instanciates a plugin from it name
    """
    def __init__(self,workdir=None, plugin_path=None):
        self.plugins = {} #name: class
        self.search_plugins()
        self._sem = Semaphore()

    def search_plugins(self, directory=None):
        """
        Search for all plugins into this directory
        """
        if not directory:
            dahu_root = os.path.dirname(os.path.abspath(__file__))
            directory = os.path.join(DAHU_ROOT, "plugins")
        if directory not in sys.path:
            sys.path.insert(0,directory)
        if os.path.isdir(directory):
            py_files = [os.path.join(directory,afile) 
                        for afile in os.listdir(directory)
                        if os.path.isfile(afile) and \
                        afile.endswith(".py")]
            for mod_file in py_files:
                try:
                    module = __import__(py_files)
                except ImportError:
                    pass
                else:
                    pass
                    
##TODO

    def __call__(self, name):
        """
        create a plugin instance from its name
        """
        if name in self.plugins:
            return self.plugins[name]()
        with self._sem:
            if "." in name:
                splitted = name.split(".")
                module_name = ".".join(splitted[:-1])
                class_name = splitted[-1]
            else:
                logger.error("plugin name have to be fully qualified")
                module_name = "dahu.plugins"
                class_name = name
            if module_name in sys.modules:
                module = sys.modules[module_name]
            else:
                module = __import__(module_name)
            assert module_name.startswith(module.__name__)
            remaining = module_name[len(module.__name__)+1:]
            if remaining:
                module =  module.__getattribute__(remaining)
            klass = module.__dict__[ class_name ]
            assert klass.IS_DAHU_PLUGIN
            self.plugins[name] = klass
        return self.plugins[name]()

        
    
    
        
plugin_factory = Factory(get_workdir()) 
