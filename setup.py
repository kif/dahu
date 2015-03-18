#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
from __future__ import with_statement, print_function

"""
Data Analysis Highly tailored for Upbl09a 
"""
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "18/03/2015"
__status__ = "development"

import os, sys, subprocess

try:
    # setuptools allows the creation of wheels
    from setuptools import setup, Command
    from setuptools.command.sdist import sdist
    from setuptools.command.build_ext import build_ext
    from setuptools.command.install_data import install_data
except ImportError:
    from distutils.core import setup, Command
    from distutils.core import Extension
    from distutils.command.sdist import sdist
    from distutils.command.build_ext import build_ext
    from distutils.command.install_data import install_data

import glob

cmdclass = {}

def get_version():
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "dahu-src"))
    import _version
    sys.path.pop(0)
    return _version.strictversion

version = get_version()

class PyTest(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.chdir("test")
        errno = subprocess.call([sys.executable, 'test_all.py'])
        if errno != 0:
            raise SystemExit(errno)
        else:
            os.chdir("..")
cmdclass['test'] = PyTest

script_files = glob.glob("scripts/*")
setup(name='dahu',
      version=version,
      author="Jérôme Kieffer (python)",
      author_email="jerome.kieffer@esrf.fr",
      description='Python lightweight, plugin based, data analysis',
      url="https://github.com/kif",
      download_url="https://github.com/kif/dahu",
      ext_package="dahu",
      scripts=script_files,
#      ext_modules=[Extension(**dico) for dico in ext_modules],
      packages=["dahu", "dahu.plugins"],
      package_dir={"dahu": "dahu-src", "dahu.plugins": "plugins" },
      test_suite="test",
      cmdclass=cmdclass,
#      data_files=data_files
      )
