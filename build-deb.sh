#!/bin/sh
export PATH=$PATH:/opt/bliss/conda/venv/dahu/bin
rm -rf deb_dist/
/usr/bin/python3 -m pip wheel . --no-cache-dir
wheel2deb default --output-dir deb_dist  --exclude numpy* 
cd deb_dist/python3-dahu*_amd64
dpkg-buildpackage -r -uc -us
cd ..
sudo dpkg -i python3-dahu*.deb
cd ..

