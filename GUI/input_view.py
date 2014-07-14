#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Project: Azimuthal integration
#             https://github.com/kif/pyFAI
#
#    Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
#
#    Principal author:       Jérôme Kieffer (Jerome.Kieffer@ESRF.eu)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""
pyFAI-calib

A tool for determining the position of a detector using a reference
sample called calibrant using Debye-Scerrer rings.

"""

__author__ = u"Jérôme Kieffer"
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "GPLv3+"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "11/07/2014"
__status__ = "development"

import logging
import os
import os.path as op
import json
from pyFAI.gui_utils import pylab, QtGui, QtCore, uic, matplotlib
from pyFAI.utils import float_, int_, str_
from pyFAI.units import hc
from pyFAI.io import is_hdf5
import pyFAI
import sys
from PyQt4.QtCore import SIGNAL
logger = logging.getLogger("calibration_view")


class InputWidget(QtGui.QWidget):
    """
    widget handling
    """
    def __init__(self, ui_file="input_widget.ui"):
        QtGui.QWidget.__init__(self)
        uic.loadUi(ui_file, self)
        self.wavelength.setValidator(QtGui.QDoubleValidator())
        self.energy.setValidator(QtGui.QDoubleValidator())
        self.pixel1.setValidator(QtGui.QDoubleValidator())
        self.pixel2.setValidator(QtGui.QDoubleValidator())

        all_detectors = pyFAI.detectors.ALL_DETECTORS.keys()
        all_detectors.sort()
        self.detector.addItems([i.capitalize() for i in all_detectors])
        self.detector.setCurrentIndex(all_detectors.index("detector"))

        self.connect(self.detector, SIGNAL("currentIndexChanged(int)"), self.detector_changed)
        self.connect(self.wavelength, SIGNAL("editingFinished ()"), self.wavelength_changed)
        self.connect(self.energy, SIGNAL("editingFinished ()"), self.energy_changed)

        # connect file selection windows
        self.connect(self.file_data, SIGNAL("clicked()"), self.select_ponifile)
        self.connect(self.file_mask_file, SIGNAL("clicked()"), self.select_maskfile)
        self.connect(self.file_dark_current, SIGNAL("clicked()"), self.select_darkcurrent)
        self.connect(self.file_flat_field, SIGNAL("clicked()"), self.select_flatfield)
        self.connect(self.file_dspacing_file, SIGNAL("clicked()"), self.select_splinefile)






    def dump(self, filename=None):
        """
        Reads all parameters and returns them as a python dictionary

        @param filename: dump the content to a file as a JSON string
        @return: dictionary

        """
        logger.debug("InputWidget.dump")
        res = {}
        res["wavelength"] = float_(self.wavelength.text()) * 1e-10

        calibrant = {"do_dspacing": bool(self.do_dspacing.isChecked())}
        calibrant["file"] = str_(self.dspacing.text()).strip()
        calibrant["calibrant"] = str_(self.detector.currentText())
        res["calibrant"] = calibrant

        detector = {"name": str_(self.detector.currentText()),
                    "file":str(self.detectorfile.text()).strip(),
                    "pixel1": float_(self.pixel1.text()),
                    "pixel2":float_(self.pixel2.text())}
        res[detector] = detector

        data = {"file": str_(self.data_file.text()).strip(),
                "hdf5": str_(self.data_file_hdf5.text()).strip() or None}
        res["data"] = data

        mask = {"file": str_(self.mask_file.text()).strip(),
                "hdf5": str_(self.mask_file_hdf5.text()).strip() or None,
                "apply":  bool(self.do_mask.isChecked())}
        res["mask"] = mask

        dark = {"file": str_(self.dark_current.text()).strip(),
                "hdf5": str_(self.dark_file_hdf5.text()).strip() or None,
                "apply":  bool(self.do_dark.isChecked())}
        res["dark"] = dark

        flat = {"file": str_(self.flat_field.text()).strip(),
                "hdf5": str_(self.flat_file_hdf5.text()).strip() or None,
                "apply":  bool(self.do_flat.isChecked())}
        res["flat"] = flat

        try:
            with open(filename, "w") as myFile:
                json.dump(res, myFile, indent=4)
        except IOError as error:
            logger.error("Error while saving config: %s" % error)
        else:
            logger.debug("Saved")

        return res

    def detector_changed(self):
        """
        Called by the UI when the combo-box value is changed
        """
        logger.debug("detector_changed")
        detector = str(self.detector.currentText()).lower()
        inst = pyFAI.detectors.detector_factory(detector)
        if inst.force_pixel:
            self.pixel1.setText(str(inst.pixel1))
            self.pixel2.setText(str(inst.pixel2))
            self.detectorfile.setText("")
        elif self.detectorfile.text():
            detectorfile = str(self.detectorfile.text()).strip()
            if op.isfile(detectorfile):
                inst.set_splineFile(detectorfile)
                self.pixel1.setText(str(inst.pixel1))
                self.pixel2.setText(str(inst.pixel2))
            else:
                logger.warning("No such spline file %s" % detectorfile)


    def energy_changed(self):
        logger.debug("energy changed")
        self.wavelength.setText(str(hc / float(self.energy.text())))

    def wavelength_changed(self):
        logger.debug("wavelength_changed")
        self.energy.setText(str(hc / float(self.wavelength.text())))

    def select_datafile(self):
        datafile = QtGui.QFileDialog.getOpenFileName()
        self.data_file.setText(datafile)
        if is_hdf5(datafile):
            pass
        #Check here for HDF5 and open next windows if needed


    def select_splinefile(self):
        logger.debug("select_splinefile")
        splinefile = str(QtGui.QFileDialog.getOpenFileName())
        if splinefile:
            try:
                self.ai.detector.set_splineFile(splinefile)
                self.pixel1.setText(str(self.ai.pixel1))
                self.pixel2.setText(str(self.ai.pixel2))
                self.splineFile.setText(self.ai.detector.splineFile or "")
            except Exception as error:
                logger.error("failed %s on %s" % (error, splinefile))

    def select_maskfile(self):
        logger.debug("select_maskfile")
        maskfile = str(QtGui.QFileDialog.getOpenFileName())
        if maskfile:
            self.mask_file.setText(maskfile or "")
            self.do_mask.setChecked(True)

    def select_darkcurrent(self):
        logger.debug("select_darkcurrent")
        darkcurrent = str(QtGui.QFileDialog.getOpenFileName())
        if darkcurrent:
            self.dark_current.setText(str_(darkcurrent))
            self.do_dark.setChecked(True)

    def select_flatfield(self):
        logger.debug("select_flatfield")
        flatfield = str(QtGui.QFileDialog.getOpenFileName())
        if flatfield:
            self.flat_field.setText(str_(flatfield))
            self.do_flat.setChecked(True)


if __name__ == "__main__":
    import fabio
    app = QtGui.QApplication(sys.argv)
    cw = InputWidget()
    cw.show()
    sys.exit(app.exec_())
