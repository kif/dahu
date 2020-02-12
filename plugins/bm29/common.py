#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

Common data structures: Sample, Ispyb
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "12/02/2020"
__status__ = "development"
version = "0.0.1"

from typing import NamedTuple
import logging
logger = logging.getLogger(__name__)


def _fromdict(cls, dico):
    "Mirror of _asdict: take the dict and populate the tuple to be returned"
    try:
        obj = cls(**dico)
    except TypeError as err:
        logger.warning("TypeError: %s", err)
        obj = cls(**{key: dico[key] for key in [i for i in cls._fields if i in dico]})
    return obj


class Sample(NamedTuple):
    """ This object represents the sample with the following representation 
      "sample": {
        "name": "bsa",
        "description": "protein description like Bovine Serum Albumin",
        "buffer": "description of buffer, pH, ...",
        "concentration": 0,
        "hplc": "column name and chromatography conditions",
        "temperature": 20,
        "temperature_env": 20},  
    """
    name: str=None
    description: str=None
    buffer: str=None
    concentration: float=None
    hplc: str=None
    storage_temp: float=None
    exposure_temp: float=None

    _fromdict = classmethod(_fromdict)
    
    
        
class Ispyb(NamedTuple):        
    url: str=None
    login: str=None
    passwd: str=None
    pyarch: str=None

    _fromdict = classmethod(_fromdict)


class EquivalentFrames(NamedTuple):
    start: int=0
    end: int=-1

def get_equivalent_frames(proba, absolute=0.1, relative=0.2):
    """This function return the start and end index of a set of equivalent data:
    Note the end-index is excluded...
    
    :param proba: 2D array with the probablility that 2 dataset are array_equivalent
    :param absolute: minimum probablity of 2 dataset to be considered equivalent
    :param absolute: minimum probablity of 2 adjacents dataset to be considered equivalent
    """
    res = []
    sizes = []
    size = len(proba)
    diag = numpy.diagonal(proba, 1) > relative
    "diag is the true if this a start point worth to be considered"
    for start, valid in enumerate(diag):
        if not valid:
            continue
        # searching for the end-point which is the first invalid
        invalid = proba[start, :] < absolute
        invalid[:start] = False  # only search for greater indices
        wend = numpy.where(invalid)[0]
        end = wend[0] if len(wend) else size
        sizes.append(end-start)
        res.append(EquivalentFrames(start, end))
    return res[numpy.argmax(sizes)]

        
            