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

from dahu.factory import optional_plugin, register

# One block per plugin: a missing dependency disables only the plugin concerned.
with optional_plugin("bm29.integratemultiframe"):
    from .integrate import IntegrateMultiframe
    register(IntegrateMultiframe, fqn="bm29.integratemultiframe")

with optional_plugin("bm29.subtractbuffer"):
    from .subtracte import SubtractBuffer
    register(SubtractBuffer, fqn="bm29.subtractbuffer")

with optional_plugin("bm29.hplc"):
    from .hplc import HPLC
    register(HPLC, fqn="bm29.hplc")

with optional_plugin("bm29.mesh"):
    from .mesh import Mesh
    register(Mesh, fqn="bm29.mesh")
