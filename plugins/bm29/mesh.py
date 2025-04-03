#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* Mesh mode: Rebuild the complete map and performs basic analysis on it.
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "03/04/2025"
__status__ = "development"
__version__ = "0.1.0"

import json
from dataclasses import dataclass

@dataclass(slots=True)
class Scan:
    """Class describing 2D mesh-scan"""
    fast_motor_name: str = "fast"
    fast_motor_start: float = 0.0
    fast_motor_stop: float = 0.0
    fast_motor_step: int = 0  # bliss-like steps, actually step+1 points
    slow_motor_name: str = "slow"
    slow_motor_start: float = 0.0
    slow_motor_stop: float = 0.0
    slow_motor_step: int = 0  # bliss-like steps, actually step+1 points
    backnforth: bool = False

    def __repr__(self):
        return json.dumps(self.as_dict(), indent=4)

    def as_dict(self):
        """Like asdict, without extra features:

        :return: dict which can be JSON-serialized 
        """
        dico = {}
        for key, value in asdict(self).items():
            dico[key] = value
        return dico

    @classmethod
    def from_dict(cls, dico):
        """Alternative constructor,
        :param dico: dict with the config
        :return: instance of the dataclass
        """
        to_init = dico
        self = cls(**to_init)
        return self

    def save(self, filename):
        """Dump the content of the dataclass as JSON file"""
        with open(filename, "w") as w:
            w.write(json.dumps(self.as_dict(), indent=2))




def input_from_master(master_file):
    """Convert a bliss masterfile containing a single NXentry with a 2D scan into a set of plugins to be launched

    :param master_file: path to a bliss masterfile
    :return: list of json dicts.
    """
    pass


class HPLC(Plugin):
    """ Rebuild the complete map and perform basic analysis on it.
    
        Typical JSON file:
    {
      "integrated_files": ["img_001.h5", "img_002.h5"],
      "output_file": "mesh.h5"
      "ispyb": {
        "url": "http://ispyb.esrf.fr:1234",
        "pyarch": "/data/pyarch/mx1234/sample", 
        "measurement_id": -1,
        "collection_id": -1
       },
       "nmf_components": 5, 
       "scan": {
            "fast_motor_name": "chipz",
            "fast_motor_start": -2.2,
            "fast_motor_stop": -3.2,
            "fast_motor_step": 3, 
            "slow_motor_name": "chipy",
            "slow_motor_start": -5.2,
            "slow_motor_stop": -12.2,
            "slow_motor_step": 7, 
            "backnforth": False,
            }
      "wait_for": [jobid_img001, jobid_img002],
      "plugin_name": "bm29.mesh"
    } 
    """
