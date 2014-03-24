#!/usr/bin/env python
# -*- coding: utf8 -*-
#
from __future__ import with_statement, print_function
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
    recompiled = False
    name = "dahu"
#    image_home = os.path.join(test_home, "testimages")
#    if not os.path.isdir(image_home):
#        os.makedirs(image_home)
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
#    opencl = os.path.join(os.path.dirname(test_home), "openCL")
#    for clf in os.listdir(opencl):
#        if clf.endswith(".cl") and clf not in os.listdir(os.path.join(dahu_home, "dahu")):
#            copy(os.path.join(opencl, clf), os.path.join(dahu_home, "dahu", clf))
    dahu = imp.load_module(*((name,) + imp.find_module(name, [dahu_home])))
    sys.modules[name] = dahu
    logger.info("dahu loaded from %s" % dahu.__file__)

    @classmethod
    def deep_reload(cls):
        logger.info("Loading dahu")
        cls.dahu = None
        dahu = None
        sys.path.insert(0, cls.dahu_home)
        for key in sys.modules.copy():
            if key.startswith("dahu"):
                sys.modules.pop(key)

        import dahu
        cls.dahu = dahu
        logger.info("dahu loaded from %s" % dahu.__file__)
        sys.modules[cls.name] = dahu
        return cls.dahu

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
#                    opencl = os.path.join(os.path.dirname(cls.test_home), "openCL")
#                    for clf in os.listdir(opencl):
#                        if clf.endswith(".cl") and clf not in os.listdir(os.path.join(cls.dahu_home, "dahu")):
#                            copy(os.path.join(opencl, clf), os.path.join(cls.dahu_home, "dahu", clf))
                    cls.dahu = cls.deep_reload()
                    cls.recompiled = True


#
#    @classmethod
#    def timeoutDuringDownload(cls, imagename=None):
#            """
#            Function called after a timeout in the download part ...
#            just raise an Exception.
#            """
#            if imagename is None:
#                imagename = "2252/testimages.tar.bz2 unzip it "
#            raise RuntimeError("Could not automatically \
#                download test images!\n \ If you are behind a firewall, \
#                please set both environment variable http_proxy and https_proxy.\
#                This even works under windows ! \n \
#                Otherwise please try to download the images manually from \n %s/%s and put it in in test/testimages." % (cls.url_base, imagename))

#
#
#    @classmethod
#    def getimage(cls, imagename):
#        """
#        Downloads the requested image from Forge.EPN-campus.eu
#        @param: name of the image.
#        For the RedMine forge, the filename contains a directory name that is removed
#        @return: full path of the locally saved file
#        """
#        baseimage = os.path.basename(imagename)
#        logger.info("UtilsTest.getimage('%s')" % baseimage)
#        fullimagename = os.path.abspath(os.path.join(cls.image_home, baseimage))
#        if not os.path.isfile(fullimagename):
#            logger.info("Trying to download image %s, timeout set to %ss"
#                          % (imagename, cls.timeout))
#            dictProxies = {}
#            if "http_proxy" in os.environ:
#                dictProxies['http'] = os.environ["http_proxy"]
#                dictProxies['https'] = os.environ["http_proxy"]
#            if "https_proxy" in os.environ:
#                dictProxies['https'] = os.environ["https_proxy"]
#            if dictProxies:
#                proxy_handler = urllib2.ProxyHandler(dictProxies)
#                opener = urllib2.build_opener(proxy_handler).open
#            else:
#                opener = urllib2.urlopen
#
##           Nota: since python2.6 there is a timeout in the urllib2
#            timer = threading.Timer(cls.timeout + 1, cls.timeoutDuringDownload, args=[imagename])
#            timer.start()
#            logger.info("wget %s/%s" % (cls.url_base, imagename))
#            if sys.version > (2, 6):
#                data = opener("%s/%s" % (cls.url_base, imagename),
#                              data=None, timeout=cls.timeout).read()
#            else:
#                data = opener("%s/%s" % (cls.url_base, imagename),
#                              data=None).read()
#            timer.cancel()
#            logger.info("Image %s successfully downloaded." % baseimage)
#
#            try:
#                open(fullimagename, "wb").write(data)
#            except IOError:
#                raise IOError("unable to write downloaded \
#                    data to disk at %s" % cls.image_home)
#
#            if not os.path.isfile(fullimagename):
#                raise RuntimeError("Could not automatically \
#                download test images %s!\n \ If you are behind a firewall, \
#                please set both environment variable http_proxy and https_proxy.\
#                This even works under windows ! \n \
#                Otherwise please try to download the images manually from \n%s/%s" % (imagename, cls.url_base, imagename))
#
#        return fullimagename


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
    mylogger.setLevel(level)
    mylogger.debug("tests loaded from file: %s" % basename)
    if force_build:
        UtilsTest.forceBuild(force_remove)
    else:
        UtilsTest.deep_reload()
    return mylogger
