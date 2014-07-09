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
__date__ = "09/07/2014"
__status__ = "development"

import fabio
import logging
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import os
import pyFAI
from pyFAI.gui_utils import pylab, QtGui, QtCore, uic
import sys


SIGNAL = pyFAI.gui_utils.QtCore.SIGNAL

logger = logging.getLogger("calibration_view")


class CalibrationWindow(QtGui.QMainWindow):
    """
    Order on the layers in the right tab
     
    mask > data > massif > solid_angle
    """
    ZORDER = {"contour":5,
              "point": 4,
              "mask": 3,
              "data": 2,
              "massif": 1,
              "solidangle": 0,
              }

    def __init__(self, imagename=None):
        QtGui.QWidget.__init__(self)
        self.imagename = imagename
        uic.loadUi("calibration_main.ui", self)
        self.connect(self.actionAbout_calibrate, SIGNAL("triggered()"), self.on_about)
        self.dpi = 100
        self.fig = self.canvas = self.mpl_toolbar = self.pix_coords_label = None
        self.data_axes = self.mask_axes = self.massif_axes = self.contour_axes = self.solid_angle_axis = None
        self.data = self.massif = self.solid_angle = None
        self.display_widget = None
        self.build_right_frame()

    def on_about(self):
        msg = [__doc__,
               "",
               "Version date: \t%s" % __date__,
               "PyFAI version: \t%s" % pyFAI.version,
               "FabIO version: \t%s" % fabio.version,
               "Author: \t\t%s" % __author__,
               "Copyright: \t\t%s" % __copyright__,
               "License: \t\t%s" % __license__]
        QtGui.QMessageBox.about(self, "About calibrate", os.linesep.join(msg))

    def build_right_frame(self):
        "build the right frame that contains matplotlib widgets"

        self.fig = Figure(dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.image_frame)
        # Create the navigation toolbar, tied to the canvas
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.image_frame, coordinates=False)

        self.data_axes = self.fig.add_subplot(111)
        self.data_axes.set_visible(False)
        # Bind the 'pick' event for clicking on one of the bars
        self.canvas.mpl_connect('motion_notify_event', self.on_pick)

        self.pix_coords_label = QtGui.QLabel("x= None , y= None , i= None ", self)
        self.mpl_toolbar.addWidget(self.pix_coords_label)
        self.display_widget = uic.loadUi("display_widget.ui")
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.mpl_toolbar, alignment=QtCore.Qt.AlignVCenter)
        vbox.addWidget(self.canvas, alignment=QtCore.Qt.AlignVCenter)
        vbox.addWidget(self.display_widget, alignment=QtCore.Qt.AlignVCenter)
        self.image_frame.setLayout(vbox)
        # few signals about those new widgets:
        self.connect(self.display_widget.show_data, SIGNAL("stateChanged(int)"), self.toggle_show_data)

    def on_pick(self, event):
        if event.inaxes and self.data is not None and self.data.any():
            if int(event.xdata) <= self.data.shape[1] and int(event.ydata) <= self.data.shape[0]:
                x = int(round(event.xdata))
                y = int(round(event.ydata))
                i = self.data[y, x]
                self.pix_coords_label.setText("x=%6d, y=%6d, I=%6g" % (x, y, i))
            else:
                self.pix_coords_label.setText("x= None , y= None , I= None ")
        else:
            self.pix_coords_label.setText("x= None , y= None , I= None ")

    def set_data(self, data, display=True):
        """
        Display an array in the  data layer  
        @param data: the numpy array with the data in it
        @param displa: shall the data  
        """
        self.data_axes.imshow(data)
        self.data_axes.set_visible(display)
        self.canvas.draw()
        if display:
            self.display_widget.show_data.setChecked(True)

    def toggle_show_data(self):
        if self.data_axes is not None:
            self.data_axes.set_visible(self.display_widget.show_data.isChecked())
            self.canvas.draw()

    def toggle_show_mask(self):
        if self.mask_axes is not None:
            self.mask_axes.set_visible(self.display_widget.show_mask.isChecked())
            self.canvas.draw()

    def toggle_show_data(self):
        if self.data_axes is not None:
            self.data_axes.set_visible(self.display_widget.show_data.isChecked())
            self.canvas.draw()


if __name__ == "__main__":
    filename = os.path.join(os.environ["HOME"], "workspace", "pyFAI", "test",
                            "testimages", "Pilatus1M.edf")
    app = QtGui.QApplication(sys.argv)
    cw = CalibrationWindow()

    cw.show()
    cw.set_data(fabio.open(filename).data)
    sys.exit(app.exec_())
