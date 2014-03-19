#!/usr/bin/env python
# -*- coding: utf8 -*-
#
from __future__ import with_statement, print_function

"""
Data Analysis Highly tailored for Upbl09a 
"""
__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20140318"
__status__ = "development"
version = "0.1"
import os, sys
from distutils.core import setup, Extension, Command
import glob

script_files = glob.glob("scripts/*")
setup(name='dahu',
      version=version,
      author="Jérôme Kieffer (python)",
      author_email="jerome.kieffer@esrf.fr",
      description='Python lightweight plugin pipelines',
      url="https://github.com/kif",
      download_url="https://github.com/kif",
      ext_package="dahu",
      scripts=script_files,
#      ext_modules=[Extension(**dico) for dico in ext_modules],
      packages=["dahu", "dahu.plugins"],
      package_dir={"dahu": "dahu-src", "dahu.plugins": "plugins" },
      test_suite="test",
#      cmdclass=cmdclass,
#      data_files=data_files
      )
