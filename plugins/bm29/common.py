#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

Common data structures: Sample, Ispyb
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "15/05/2020"
__status__ = "development"
version = "0.0.1"

import os
from collections import namedtuple
from typing import NamedTuple
import json
import logging
logger = logging.getLogger("bm29.common")
import numpy
from dahu.cache import DataCache
from hdf5plugin import Bitshuffle, Zfp
import pyFAI, pyFAI.units
from pyFAI.method_registry import IntegrationMethod
import fabio
if pyFAI.version_info < (0, 20):
    from .nexus import Nexus, get_isotime
else:
    from pyFAI.io import Nexus, get_isotime
#cmp contains the compression options, shared by all plugins. Used mainly for images 
cmp_int = Bitshuffle()
cmp_float = Zfp(rate=8)
cmp = cmp_float

# Constants associated to the azimuthal integrator to be used in all plugins:
polarization_factor = 0.9
method = IntegrationMethod.select_method(1, "no", "csr", "opencl")[0]

#This cache contains azimuthal integrators shared between the different plugins
KeyCache = namedtuple("KeyCache", "npt unit poni mask energy")
shared_cache = DataCache(10)

def get_integrator(keycache):
    "retrieve or build an azimuthal integrator based on the keycache provided "
    #key = KeyCache(self.npt, self.unit, self.poni, self.mask, self.energy)
    if keycache in shared_cache:
        ai = shared_cache[keycache]
    else:
        if os.path.exists(keycache.poni):
            ai = pyFAI.load(keycache.poni)
        else:
            logger.warning("Poni file does not exist, try to consider it as a JSON-dict")
            config = json.loads(keycache.poni)
            ai = pyFAI.load(config)
            
        ai.wavelength = 1e-10 * pyFAI.units.hc / keycache.energy
        if keycache.mask:
            mask = numpy.logical_or(fabio.open(keycache.mask).data, ai.detector.mask).astype("int8")
            ai.detector.mask = mask
        shared_cache[keycache] = ai
    return ai

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
    temperature_env: float=None
    temperature: float=None

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

    Note: the end-index is excluded, as usual in Python
    
    :param proba: 2D array with the probablility that 2 dataset are array_equivalent
    :param absolute: minimum probablity of 2 dataset to be considered equivalent
    :param absolute: minimum probablity of 2 adjacents dataset to be considered equivalent
    :return: 2-tuple with start and end for the region of equivalent data
    """
    res = []
    sizes = []
    size = len(proba)
    ext_diag = numpy.zeros(size+1, dtype=numpy.int16)
    delta = numpy.zeros(size+1, dtype=numpy.int16)
    ext_diag[0] = 1
    ext_diag[1:-1] = numpy.diagonal(proba, 1) >= relative
    delta[0] = ext_diag[1]
    delta[1:] = ext_diag[1:] - ext_diag[:-1]
    start = numpy.where(delta>0)[0]
    end = numpy.where(delta<0)[0]
    for start_i, end_i in zip(start, end):
        for start_j in range(start_i, end_i):
            # searching for the end-point which is the first invalid
            invalid = proba[start_j, :] < absolute
            invalid[:start_j] = False  # only search for greater indices
            wend = numpy.where(invalid)[0]
            end_j = min(end_i, wend[0] if len(wend) else size)
            sizes.append(end_j - start_j)
            res.append((start_j, end_j))
    return res[numpy.argmax(sizes)]
