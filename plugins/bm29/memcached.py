"""Data Analysis plugin for BM29: BioSaxs

Connection to  Memcached
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "17/09/2026"
__status__ = "development"
__version__ = "0.3.0"

import socket

try:
    import memcache
except (ImportError, ModuleNotFoundError):
    memcache = None

SERVER = "localhost"


def to_memcached(dico):
    rc = {}
    if memcache is not None:
        mc = memcache.Client([(SERVER, 11211)])
        rc["server"] = socket.getfqdn()+":11211"
        for k, v in dico.items():
            if len(k)>250:
                k = k[-250:]
            rc[k] = mc.set(k, v)
    return rc
