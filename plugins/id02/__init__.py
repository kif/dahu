"""Set of plugins for ID02/Trusaxs

List of plugins exposed:
* Metadata saving (C216):         id02.Metadata
* single detector processing:     id02.SingleDetector
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "16/03/2020"
__status__ = "development"
__version__ = "1.0.0"

from dahu.factory import register
from .metadata import Metadata
from .single_detector import SingleDetector

register(Metadata, fqn="id02.metadata")
register(SingleDetector, fqn="id02.singledetector")
