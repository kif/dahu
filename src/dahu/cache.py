#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""
Data Analysis RPC server over Tango: 

Class Cache for storing the data in a Borg  
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "22/02/2022" 
__status__ = "production"

import os
import logging
from threading import Semaphore
logger = logging.getLogger("dahu.cache")


class DataCache(dict):
    """
    Class behaves like a dict with a finite size, used to cache some data for some time 
    then discard them when other things are stored.
    
    This class can be configured as a borg (singleton-like): 
    It always returns the same values regardless to the instance of the object
    """
    __shared_state = {}
    __data_initialized = False

    def __init__(self, max_size=10, borg=True):
        """
        Constructor of DataCache
        :param max_size: number of element to keep in memory
        :param borg: set to false to have the behavour of a normal class. By default, this is a Borg, all instances have the same content. 
        """
        if borg:
            self.__dict__ = self.__shared_state
            if DataCache.__data_initialized is False:
                DataCache.__data_initialized = True
                logger.debug("DataCache.__init__: initalization of the Borg")
                self.ordered = []
                self.dict = {}
                self.max_size = max_size
                self._sem = Semaphore()
            else:
                logger.debug("DataCache.__init__: Borg already initialized")
        else:
            logger.debug("DataCache.__init__: initalization of a standard class")
            self.ordered = []
            self.dict = {}
            self.max_size = max_size
            self._sem = Semaphore()
            
    def __repr__(self):
        """
        """
        out = ["{"]
        for key in self.ordered:
            out.append(" '%s': %s," % (key, self.dict[key]))
        out.append("}")
        return os.linesep.join(out)

    def __setitem__(self, key, value):
        """
        x.__setitem__(i, y) <==> x[i]=y
        """
        with self._sem:
            logger.debug("DataCache.__setitem__: %s", key)
            self.dict[key] = value
            if key in self.ordered:
                index = self.ordered.index(key)
                self.ordered.pop(index)
            if len(self.ordered) > self.max_size:
                firstKey = self.ordered.pop(0)
                logger.debug("Removing from cache: %s", firstKey)
                self.dict.pop(firstKey)
            self.ordered.append(key)

    def __getitem__(self, key):
        """
        x.__getitem__(y) <==> x[y]
        """
        with self._sem:
            logger.debug("DataCache.__setitem__: %s", key)
            index = self.ordered.index(key)
            self.ordered.pop(index)
            self.ordered.append(key)
            return self.dict[key]

    def __contains__(self, key):
        """
        D.__contains__(k) -> True if D has a key k, else False
        """
        return key in self.dict
    has_key = __contains__

    def __len__(self):
        """
        Returns the length of the object
        """
        return len(self.ordered)

    def get(self, key, default=None):
        """
        get method with default answer implemented
        """
        if key in self.ordered:
            return self.__getitem__(key)
        elif default is not None:
            self.__setitem__(key, default)
            return default

    def keys(self):
        """
        Returns the list of keys, ordered
        """
        logger.debug("DataCache.keys")
        return self.ordered[:]

    def pop(self, key):
        """
        Remove a key for the dictionary and return it's value
        """
        with self._sem:
            logger.debug("DataCache.pop %s", key)
            try:
                index = self.ordered.index(key)
            except:
                raise KeyError
            self.ordered.pop(index)
            myData = self.dict.pop(key)
        return myData
