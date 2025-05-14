"""Set of plugins for BM29/BioSaxs

List of plugins exposed:

* bm29.IntegrateMultiframe
* bm29.SubtractBuffer
* bm29.hplc
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "05/05/2025"
__status__ = "development"
__version__ = "0.2.0"

from dahu.factory import register
from .integrate import IntegrateMultiframe
from .subtracte import SubtractBuffer
from .hplc import HPLC
from .mesh import Mesh
register(IntegrateMultiframe, fqn="bm29.integratemultiframe")
register(SubtractBuffer, fqn="bm29.subtractbuffer")
register(HPLC, fqn="bm29.hplc")
register(Mesh, fqn="bm29.mesh")