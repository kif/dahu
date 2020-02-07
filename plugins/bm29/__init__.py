"""Set of plugins for BM29/BioSaxs

List of plugins exposed:

* bm29.IntegrateMultiframe
"""
from dahu.factory import register
from .integrate import IntegrateMultiframe
register(IntegrateMultiframe, fqn="bm29.integratemultiframe")