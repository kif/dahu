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
__date__ = "20140611"
__status__ = "development"

import os, sys
from distutils.core import setup, Extension, Command
import glob

version = [eval(l.split("=")[1]) for l in open(os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "dahu-src", "__init__.py"))
    if l.strip().startswith("version")][0]

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
#      cmdclass=cmdclass,
#      data_files=data_files
      )
