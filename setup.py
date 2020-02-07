#!/usr/bin/env python
# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2013-2016 European Synchrotron Radiation Facility
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/

from __future__ import with_statement, print_function

"""
Data Analysis Highly tailored for Upbl09a 
"""
__authors__ = ["Jérôme Kieffer", "Thomas Vincent"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "06/02/2020"
__status__ = "development"


import sys
import os
import platform
import glob

from numpy.distutils.misc_util import Configuration

try:
    # setuptools allows the creation of wheels
    from setuptools import setup, Command
    from setuptools.command.build_py import build_py as _build_py
    from setuptools.command.sdist import sdist
    from setuptools.command.build_ext import build_ext
    from setuptools.command.install_data import install_data
    from setuptools.command.sdist import sdist
except ImportError:
    from numpy.distutils.core import setup, Command
    from distutils.command.build_py import build_py as _build_py
    from distutils.core import Extension
    from distutils.command.sdist import sdist
    from distutils.command.build_ext import build_ext
    from distutils.command.install_data import install_data
    from distutils.command.sdist import sdist

PROJECT = "dahu"
cmdclass = {}

def get_version():
    import version
    return version.strictversion



def get_readme():
    dirname = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(dirname, "README.rst"), "r") as fp:
        long_description = fp.read()
    return long_description


classifiers = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Environment :: MacOS X",
               "Environment :: Win32 (MS Windows)",
               "Environment :: X11 Applications :: Qt",
               "Intended Audience :: Education",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)",
               "Natural Language :: English",
               "Operating System :: MacOS",
               "Operating System :: Microsoft :: Windows",
               "Operating System :: POSIX",
               "Programming Language :: Cython",
               "Programming Language :: Python :: 2.7",
               "Programming Language :: Python :: 3.4",
               "Programming Language :: Python :: 3.5",
               "Programming Language :: Python :: Implementation :: CPython",
               "Topic :: Scientific/Engineering :: Physics",
               "Topic :: Software Development :: Libraries :: Python Modules",
               ]


# ########## #
# version.py #
# ########## #

class build_py(_build_py):
    """
    Enhanced build_py which copies version.py to <PROJECT>._version.py
    """
    def find_package_modules(self, package, package_dir):
        modules = _build_py.find_package_modules(self, package, package_dir)
        if package == PROJECT:
            modules.append((PROJECT, '_version', 'version.py'))
        return modules


cmdclass['build_py'] = build_py


########
# Test #
########

class PyTest(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess
        errno = subprocess.call([sys.executable, 'run_tests.py', '-i'])
        if errno != 0:
            raise SystemExit(errno)
cmdclass['test'] = PyTest

# ################### #
# build_doc commandes #
# ################### #

try:
    import sphinx
    import sphinx.util.console
    sphinx.util.console.color_terminal = lambda: False
    from sphinx.setup_command import BuildDoc
except ImportError:
    sphinx = None
else:
    # i.e. if sphinx:
    class build_doc(BuildDoc):

        def run(self):
            # make sure the python path is pointing to the newly built
            # code so that the documentation is built on this and not a
            # previously installed version

            build = self.get_finalized_command('build')
            sys.path.insert(0, os.path.abspath(build.build_lib))

#             # Copy .ui files to the path:
#             dst = os.path.join(
#                 os.path.abspath(build.build_lib), "silx", "gui")
#             if not os.path.isdir(dst):
#                 os.makedirs(dst)
#             for i in os.listdir("gui"):
#                 if i.endswith(".ui"):
#                     src = os.path.join("gui", i)
#                     idst = os.path.join(dst, i)
#                     if not os.path.exists(idst):
#                         shutil.copy(src, idst)

            # Build the Users Guide in HTML and TeX format
            for builder in ('html', 'latex'):
                self.builder = builder
                self.builder_target_dir = os.path.join(self.build_dir, builder)
                self.mkpath(self.builder_target_dir)
                BuildDoc.run(self)
            sys.path.pop(0)
    cmdclass['build_doc'] = build_doc

# ############################# #
# numpy.distutils Configuration #
# ############################# #

def configuration(parent_package='', top_path=None):
    """Recursive construction of package info to be used in setup().

    See http://docs.scipy.org/doc/numpy/reference/distutils.html#numpy.distutils.misc_util.Configuration
    """  # noqa
    config = Configuration(None, parent_package, top_path)
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
        quiet=True)
    config.add_subpackage(PROJECT)
    return config


config = configuration()

################################################################################
# Debian source tree
################################################################################


class sdist_debian(sdist):
    """
    Tailor made sdist for debian
    * remove auto-generated doc
    * remove cython generated .c files
    """
    def prune_file_list(self):
        sdist.prune_file_list(self)
        to_remove = ["doc/build", "doc/pdf", "doc/html", "pylint", "epydoc"]
        print("Removing files for debian")
        for rm in to_remove:
            self.filelist.exclude_pattern(pattern="*", anchor=False, prefix=rm)

        # this is for Cython files specifically: remove C & html files
        search_root = os.path.dirname(os.path.abspath(__file__))
        for root, _, files in os.walk(search_root):
            for afile in files:
                if os.path.splitext(afile)[1].lower() == ".pyx":
                    base_file = os.path.join(root, afile)[len(search_root) + 1:-4]
                    self.filelist.exclude_pattern(pattern=base_file + ".c")
                    self.filelist.exclude_pattern(pattern=base_file + ".cpp")
                    self.filelist.exclude_pattern(pattern=base_file + ".html")

    def make_distribution(self):
        self.prune_file_list()
        sdist.make_distribution(self)
        dest = self.archive_files[0]
        dirname, basename = os.path.split(dest)
        base, ext = os.path.splitext(basename)
        while ext in [".zip", ".tar", ".bz2", ".gz", ".Z", ".lz", ".orig"]:
            base, ext = os.path.splitext(base)
        if ext:
            dest = "".join((base, ext))
        else:
            dest = base
        sp = dest.split("-")
        base = sp[:-1]
        nr = sp[-1]
        debian_arch = os.path.join(dirname, "-".join(base) + "_" + nr + ".orig.tar.gz")
        os.rename(self.archive_files[0], debian_arch)
        self.archive_files = [debian_arch]
        print("Building debian .orig.tar.gz in %s" % self.archive_files[0])

cmdclass['debian_src'] = sdist_debian


# ##### #
# setup #
# ##### #

setup_kwargs = config.todict()
script_files = glob.glob("scripts/*")
install_requires = ["numpy"]
setup_requires = ["numpy"]

setup_kwargs.update(name=PROJECT,
                    version=get_version(),
                    url="https://github.com/kif/dahu",
                    author="Jérôme Kieffer",
                    author_email="silx@esrf.fr",
                    classifiers=classifiers,
                    description="Python lightweight, plugin based, data analysis",
                    long_description=get_readme(),
                    install_requires=install_requires,
                    setup_requires=setup_requires,
                    cmdclass=cmdclass,
                    package_data={'silx.resources': [
                        # Add here all resources files
                        'gui/icons/*.png',
                        ]},
                    zip_safe=False,
                    packages=["dahu", "dahu.plugins", "dahu.test"],
                    package_dir={"dahu": "dahu",
                                 "dahu.plugins": "plugins",
                                 "dahu.test": "dahu/test"},
                    test_suite="test",
                    scripts=script_files,
                    )

setup(**setup_kwargs)
