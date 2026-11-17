"""Set of plugins for ID02/Trusaxs

List of plugins exposed:
* Metadata saving (C216):         id02.Metadata
* single detector processing:     id02.SingleDetector
* XPCS pixel correlation          id02.xpcs
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "25/03/2020"
__status__ = "development"
__version__ = "1.0.0"

from dahu.factory import optional_plugin, register

# One block per plugin: a missing dependency disables only the plugin concerned.
with optional_plugin("id02.metadata"):
    from .metadata import Metadata
    register(Metadata, fqn="id02.metadata")

with optional_plugin("id02.singledetector"):
    from .single_detector import SingleDetector
    register(SingleDetector, fqn="id02.singledetector")

with optional_plugin("id02.xpcs"):
    from .xpcs import XPCS
    register(XPCS, fqn="id02.xpcs")
