#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

Common data structures: Sample, Ispyb
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "10/02/2020"
__status__ = "development"
version = "0.0.1"

from typing import NamedTuple
import logging
logger=logging.getLogger(__name__)


class Sample(NamedTuple):
    """ This object represents the sample with the following representation 
    {
        "code": "bsa",
        "comments": "protein name",
        "buffer": "description of buffer, pH, ...",
        "concentration": 0,
        "hplc": "column name and chromatography conditions",
        "storage_temp": 20,
        "exposure_temp": 20},
    """
    code: str=None
    comments: str=None
    buffer: str=None
    concentration: float=None
    hplc: str=None
    storage_temp: float=None
    exposure_temp: float=None

    @classmethod
    def _fromdict(cls, dico):
        "Mirror of _asdict: take the dict and populate the tuple to be returned"
        try:
            obj = cls(**dico)
        except TypeError as err:
            logger.warning("TypeError: %s", err)
            intersection = [i for i in cls._fields if i in dico]
            obj = cls(**{key: dico[key] for key in intersection})
        return obj
    
        
class Ispyb(NamedTuple):        
    url: str=None
    login: str=None
    passwd: str=None
    pyarch: str=None

    @classmethod
    def _fromdict(cls, dico):
        "Mirror of _asdict: take the dict and populate the tuple to be returned"
        try:
            obj = cls(**dico)
        except TypeError as err:
            logger.warning("TypeError: %s", err)
            intersection = [i for i in cls._fields if i in dico]
            obj = cls(**{key: dico[key] for key in intersection})
        return obj
