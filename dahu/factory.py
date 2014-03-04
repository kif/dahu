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
import os
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
        if os.path.isdir(directory):
            py_files = [os.path.join(directory,afile) 
                        for afile in os.listdir(directory)
                        if os.path.isfile(afile) and \
                        afile.endswith(".py")]
#            if a
#TODO

    def __call__(self, name):
        """
        create a plugin instance from its name
        """
        
        module = my_import(_strPluginName, strModuleLocation)
        klass = module.__dict__[ name ]
        with self._sem:
            
        self.plugins
        return plugin

        
    
    
        
plugin_factory = Factory(get_workdir()) 
