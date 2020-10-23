#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

Common data structures: Sample, Ispyb
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "23/10/2020"
__status__ = "development"
version = "0.0.2"

import os
from pathlib import Path
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
from .nexus import Nexus, get_isotime
# else:
#     from pyFAI.io import Nexus, get_isotime
    
#cmp contains the compression options, shared by all plugins. Used mainly for images 
cmp = cmp_int = Bitshuffle()
cmp_float = Zfp(reversible=True) 


#This is used for NXdata plot style
SAXS_STYLE = json.dumps({"signal_scale_type": "log"},
                        indent=2, 
                        separators=(",\r\n", ":\t"))
NORMAL_STYLE = json.dumps({"signal_scale_type": "linear"},
                          indent=2, 
                          separators=(",\r\n", ":\t"))


# Constants associated to the azimuthal integrator to be used in all plugins:
polarization_factor = 0.9
method = IntegrationMethod.select_method(1, "no", "csr", "opencl")[0]

#This cache contains azimuthal integrators shared between the different plugins
KeyCache = namedtuple("KeyCache", "npt unit poni mask energy")
shared_cache = DataCache(10)

#Try to load the default login and passwd for Ispyb:
_default_passwd = {}
_ispyb_passwd_file = os.path.join(str(Path.home()), ".ispyb")
if os.path.exists(_ispyb_passwd_file):
    with open(_ispyb_passwd_file) as f:
        _default_passwd = json.load(f)


def get_integrator(keycache):
    "retrieve or build an azimuthal integrator based on the keycache provided "
    #key = KeyCache(self.npt, self.unit, self.poni, self.mask, self.energy)
    if keycache in shared_cache:
        ai = shared_cache[keycache]
    else:
        ai = pyFAI.load(keycache.poni)
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
    name: str="Unknown sample"
    description: str=None
    buffer: str=None
    concentration: float=None
    hplc: str=None
    temperature_env: float=None
    temperature: float=None

    _fromdict = classmethod(_fromdict)

    def __repr__(self):
        return f"{self.name}, {self.concentration} mg/mL in {self.buffer}"


class Ispyb(NamedTuple):        
    url: str=None
    login: str=_default_passwd.get("username")
    passwd: str=_default_passwd.get("password")
    pyarch: str=""
    collection_id: int = -1 # This is now deprecated
    measurement_id: int=-1
    run_number: int - 1

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
    ext_diag = numpy.zeros(size + 1, dtype=numpy.int16)
    delta = numpy.zeros(size + 1, dtype=numpy.int16)
    ext_diag[0] = 1
    ext_diag[1:-1] = numpy.diagonal(proba, 1) >= relative
    delta[0] = ext_diag[1]
    delta[1:] = ext_diag[1:] - ext_diag[:-1]
    start = numpy.where(delta > 0)[0]
    end = numpy.where(delta < 0)[0]
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


def create_nexus_sample(nxs, entry, sample):
    "Create a NXsample inside the NXentry"
    sample_grp = nxs.new_class(entry, sample.name, "NXsample")
    if sample.description is not None:
        sample_grp["description"] = sample.description
    if sample.concentration is not None:
        concentration_ds = sample_grp.create_dataset("concentration", data=sample.concentration)
        concentration_ds.attrs["units"] = "mg/mL"
    if sample.buffer is not None:
        buffer_ds = sample_grp.create_dataset("buffer", data=sample.buffer)
        buffer_ds.attrs["comment"] = "Buffer description"
    if sample.hplc:
        hplc_ds = sample_grp.create_dataset("hplc", data=sample.hplc)
        hplc_ds.attrs["comment"] = "Conditions for HPLC experiment"
    if sample.temperature is not None:
        tempe_ds = sample_grp.create_dataset("temperature", data=sample.temperature)
        tempe_ds.attrs["units"] = "°C"
        tempe_ds.attrs["comment"] = "Exposure temperature"
    if sample.temperature_env is not None:
        tempv_ds = sample_grp.create_dataset("temperature_env", data=sample.temperature_env)
        tempv_ds.attrs["units"] = "°C"
        tempv_ds.attrs["comment"] = "Storage temperature"
