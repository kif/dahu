#!/bin/sh
rm -rf deb_dist
export PYBUILD_DISABLE_python2=test
export PYBUILD_DISABLE_python3=test
python3 setup.py --command-packages=stdeb.command  bdist_deb
sudo dpkg -i deb_dist/python3-*.deb
