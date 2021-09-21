"""Set of plugins for BM29/BioSaxs

List of plugins exposed:

* bm29.IntegrateMultiframe
* bm29.SubtractBuffer
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "12/10/2020"
__status__ = "development"
__version__ = "0.1.0"

from dahu.factory import register
from .integrate import IntegrateMultiframe
from .subtracte import SubtractBuffer
from .hplc import HPLC
register(IntegrateMultiframe, fqn="bm29.integratemultiframe")
register(SubtractBuffer, fqn="bm29.subtractbuffer")
register(HPLC, fqn="bm29.hplc")
