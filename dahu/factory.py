#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Data Analysis Highly tailored for Upbl09a 
"""
from __future__ import with_statement, print_function, absolute_import, division

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "07/02/2020"
__status__ = "production"

import os
import sys
import os.path as op
import logging
from collections import OrderedDict
logger = logging.getLogger("dahu.factory")
from threading import Semaphore
from .utils import get_workdir, fully_qualified_name

import importlib, importlib.util

def load_source(module_name, file_path):
    "Plugin loader which does not pollute sys.module"
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    #module = importlib.util.module_from_spec(spec)
    #spec.loader.exec_module(module)
    
    module = spec.loader.load_module(spec.name)
    print(module)
    #Option: remove from sys.modules ...
    return module


dahu_root = os.path.dirname(os.path.abspath(__file__))


class Factory(object):
    """
    This is a factory, it instanciates a plugin from it name
    """
    registry = {}
    modules = {}
    plugin_dirs = OrderedDict()  # key: directory name, value=list of modules
    reg_sem = Semaphore()

    def __init__(self, workdir=None, plugin_path=None):
        """
        @param workdir: place were we are allowed to write
        @param plugin_path: places where plugins are ... in addition to the content of DAHU_PATH"
        """
        self._sem = Semaphore()
        self.workdir = workdir or "."
        # First the ones from environment
        if "DAHU_PLUGINS" in os.environ:
            for directory in os.environ["DAHU_PLUGINS"].split(os.pathsep):
                self.add_directory(directory)
        # Then the ones from configuration
        for directory in (plugin_path or []):
            self.add_directory(directory)
        # Finally the default one
        self.add_directory(os.path.join(dahu_root, "plugins"))

    def add_directory(self, directory):
        abs_dir = os.path.abspath(directory)
        if not os.path.isdir(directory):
            logger.warning("No such directory: %s" % directory)
            return
        python_files = []
        for i in os.listdir(abs_dir):
            j = op.join(abs_dir, i)
            print(i, j, op.isdir(j), op.exists(op.join(j, "__init__.py")))
            if op.isfile(j) and i.endswith(".py"):
                python_files.append(i[:-3])
            if op.isdir(j) and op.exists(op.join(j, "__init__.py")):
                python_files.append(i)
                                         
        logger.warning(" ".join(python_files))
        with self._sem:
            self.plugin_dirs[abs_dir] = python_files

    def search_plugin(self, plugin_name):
        """
        Search for a given plugins ...
        starting from the FQN package.class, 
        """
        if "." not in plugin_name:
            logger.error("plugin name have to be fully qualified, here: %s" % plugin_name)
            return
        splitted = plugin_name.split(".")
        module_name = ".".join(splitted[:-1])

        for dirname, modules in self.plugin_dirs.items():
            if module_name in modules and module_name not in self.modules:
                dst = op.join(dirname, module_name)
                if op.isdir(dst):
                    fname = op.join(dst, "__init__.py")
                elif op.isfile(dst+".py"):
                    fname = dst+".py"
                else:
                    raise RuntimeError("Unable to find module source for %s in %s"%(module_name, dirname))
                logger.info("load %s from %s",module_name, fname)
                mod = load_source(module_name, os.path.join(dirname, fname))
                with self.reg_sem:
                    self.modules[module_name] = mod

    def __call__(self, plugin_name):
        """
        create a plugin instance from its name
        
        @param plugin_name: name of the plugin as a string
        @return: plugin instance
        """
        plugin_name = plugin_name.lower()
        if plugin_name in self.registry:
            return self.registry[plugin_name]()
        with self._sem:
            self.search_plugin(plugin_name)
        if plugin_name not in self.registry:
            logger.error("Plugin directories have been searched but plugin"
                         " %s was not found" % plugin_name)
        else:
            return self.registry[plugin_name]()

    @classmethod
    def register(cls, klass, fqn=None):
        """
        Register a class as a plugin which can be instanciated.
        
        This can be used as a decorator
        
        @plugin_factor.register 
        
        @param klass: class to be registered as a plugin
        @param fqn: fully qualified name 
        @return klass
        """
        if fqn is None:
            fqn = fully_qualified_name(klass)
        logger.debug("Registering plugin %s as %s" % (klass, fqn))
        with cls.reg_sem:
            cls.registry[fqn] = klass
        return klass

plugin_factory = Factory(get_workdir())
register = plugin_factory.register
