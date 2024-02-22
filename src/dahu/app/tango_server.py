#!/usr/bin/env python
# coding: utf-8

#Example on how to launch it:
#PYTHONPATH=build/lib.linux-x86_64-2.6/ TANGO_HOST=saxs1:20000 python scripts/dahu_server.py GPU1
#Example for Jive:
#"example.square",{\"x\":5}
"""

Data analysis Tango device server ... for UPBL09a

"""
__author__ = "Jérôme Kieffer"
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "22/02/2024"
__status__ = "beta"
__docformat__ = 'restructuredtext'

import sys
import os
import tempfile
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("dahu_server")
try:
    from rfoo.utils import rconsole
    rconsole.spawn_server()
except ImportError:
    logger.debug("No socket opened for debugging. Please install rfoo")


import dahu.utils
from dahu.server import DahuDSClass, DahuDS

# set loglevel at least at INFO
if logger.getEffectiveLevel() > logging.INFO:
    logger.setLevel(logging.INFO)
from argparse import ArgumentParser
import PyTango

def main(argv=None):
    if argv is None:
        argv = sys.argv
    logger.info("Starting Dahu Tango Device Server")
    description = """Data Analysis Tango device server 
"""
    epilog = """ Provided by the Data analysis unit - ESRF 
        """
    usage = "dahu_server [-d]  tango-options"
    parser = ArgumentParser(description=description, epilog=epilog, add_help=True)
    parser.add_argument("-V", "--version", action='version', version='%(prog)s 0.0')
    parser.add_argument("-d", "--debug", dest="debug", default=False,
                        action="store_true", help="Switch to debug mode ",)
    parser.add_argument("-v", "--verbose", dest="tango_verbose", default=None,
                        help="tango trace level")
    parser.add_argument("-f", "--file", dest="tango_file", default=None,
                        help="tango log filename")
    parser.add_argument("-l", "--log", dest="dahu_log", 
                        default=os.path.join(os.environ.get("HOME",tempfile.gettempdir()),"log"),
                        help="directory where dahu stores logs ")
#    parser.add_argument("-n", "--nbcpu", dest="nbcpu", type=int,
#                  help="Maximum bumber of processing threads to be started", default=None)
    parser.add_argument(dest="tango", nargs="*", help="Tango device server options")
    options = parser.parse_args()
    dahu.utils.get_workdir(options.dahu_log)
    tangoParam = ["DahuDS"] + options.tango
    if options.tango_verbose:
        tangoParam.append(f"-v{options.tango_verbose}")
    if options.tango_file:
        tangoParam.append(f"-file={options.tango_file}")

    # Analyse arguments and options
    if options.debug:
        logger.debug("Switch logger to debug level")
        logger.setLevel(logging.DEBUG)
        logging.root.setLevel(logging.DEBUG)

    try:
        print(tangoParam)
        py = PyTango.Util(tangoParam)
        py.add_TgClass(DahuDSClass, DahuDS, 'DahuDS')
        U = py.instance() #PyTango.Util.instance()
        U.server_init()
        U.server_run()
    except PyTango.DevFailed as err:
        logger.error('PyTango --> Received a DevFailed exception: %s' % err)
        return -1
    except Exception as err:
        logger.error('PyTango --> An unforeseen exception occurred....%s' % err)
        return -1
    return 0


if __name__ == '__main__':
    sys.exit(main())
