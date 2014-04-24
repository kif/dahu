#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
from __future__ import with_statement, print_function, division
__author__ = "Jérôme Kieffer"
__contact__ = "jerome.kieffer@esrf.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20140324"

import os, imp, sys, subprocess, threading
import distutils.util
import logging
import urllib2
import bz2
import gzip
import numpy
import shutil
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger("utilstest")


from argparse import ArgumentParser

parser = ArgumentParser(usage="Tests for dahu")
parser.add_argument("-d", "--debug", dest="debug", help="run in debugging mode",
                  default=False, action="store_true")
parser.add_argument("-i", "--info", dest="info", help="run in more verbose mode ",
                  default=False, action="store_true")
parser.add_argument("-f", "--force", dest="force", help="force the build of the library",
                  default=False, action="store_true")
parser.add_argument("-r", "--really-force", dest="remove",
                  help="remove existing build and force the build of the library",
                  default=False, action="store_true")
options = parser.parse_args()

def copy(infile, outfile):
    "link or copy file according to the OS"
    if "symlink" in dir(os):
        os.symlink(infile, outfile)
    else:
        shutil.copy(infile, outfile)

class UtilsTest(object):
    """
    Static class providing useful stuff for preparing tests.
    """
    timeout = 60        #timeout in seconds for downloading images
    url_base = "http://upload.wikimedia.org"
    #Nota https crashes with error 501 under windows.
#    url_base = "https://forge.epn-campus.eu/attachments/download"
    test_home = os.path.dirname(os.path.abspath(__file__))
    sem = threading.Semaphore()
    reloaded = False
    recompiled = False
    name = "dahu"
    platform = distutils.util.get_platform()
    architecture = "lib.%s-%i.%i" % (platform,
                                    sys.version_info[0], sys.version_info[1])
    dahu_home = os.path.join(os.path.dirname(test_home),
                                        "build", architecture)
    logger.info("dahu Home is: " + dahu_home)
    if "dahu" in sys.modules:
        logger.info("dahu module was already loaded from  %s" % sys.modules["dahu"])
        dahu = None
        sys.modules.pop("dahu")

    if not os.path.isdir(dahu_home):
        with sem:
            if not os.path.isdir(dahu_home):
                logger.warning("Building dahu to %s" % dahu_home)
                p = subprocess.Popen([sys.executable, "setup.py", "build"],
                                 shell=False, cwd=os.path.dirname(test_home))
                logger.info("subprocess ended with rc= %s" % p.wait())
                recompiled = True

    @classmethod
    def deep_reload(cls):
        logger.info("Loading dahu")
        cls.dahu = None
        dahu = None
        if sys.path[0] != cls.dahu_home:
            sys.path.insert(0, cls.dahu_home)
        if "dahu" in sys.modules:
            logger.info("dahu module was already loaded from  %s" % sys.modules["dahu"])
            cls.dahu = None
            sys.modules.pop("dahu")

        for key in sys.modules.copy():
            if key.startswith("dahu"):
                sys.modules.pop(key)

        dahu = imp.load_module(*((cls.name,) + imp.find_module(cls.name, [cls.dahu_home])))
        sys.modules[cls.name] = dahu
        logger.info("dahu loaded from %s" % dahu.__file__)
        cls.dahu = dahu
        cls.reloaded = True
        return dahu

    @classmethod
    def forceBuild(cls, remove_first=True):
        """
        force the recompilation of dahu
        """
        if not cls.recompiled:
            with cls.sem:
                if not cls.recompiled:
                    logger.info("Building dahu to %s" % cls.dahu_home)
                    if "dahu" in sys.modules:
                        logger.info("dahu module was already loaded from  %s" % sys.modules["dahu"])
                        cls.dahu = None
                        sys.modules.pop("dahu")
                    if remove_first:
                        recursive_delete(cls.dahu_home)
                    p = subprocess.Popen([sys.executable, "setup.py", "build"],
                                     shell=False, cwd=os.path.dirname(cls.test_home))
                    logger.info("subprocess ended with rc= %s" % p.wait())
                    cls.deep_reload()
                    cls.recompiled = True




def recursive_delete(strDirname):
    """
    Delete everything reachable from the directory named in "top",
    assuming there are no symbolic links.
    CAUTION:  This is dangerous!  For example, if top == '/', it
    could delete all your disk files.
    @param strDirname: top directory to delete
    @type strDirname: string
    """
    for root, dirs, files in os.walk(strDirname, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(strDirname)

def getLogger(filename=__file__):
    """
    small helper function that initialized the logger and returns it
    """
    basename = os.path.basename(os.path.abspath(filename))
    basename = os.path.splitext(basename)[0]
    force_build = False
    force_remove = False
    level = logging.WARN
    if options.debug:
        level = logging.DEBUG
    elif options.info:
        level = logging.INFO
    if options.force:
        force_build = True
    if options.remove:
        force_remove = True
        force_build = True
    mylogger = logging.getLogger(basename)
    logger.setLevel(level)
    mylogger.setLevel(level)
    mylogger.debug("tests loaded from file: %s" % basename)
    if force_build:
        UtilsTest.forceBuild(force_remove)
    return mylogger
