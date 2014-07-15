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
__date__ = "10/07/2014"
__status__ = "development"

import fabio
import logging
import matplotlib
matplotlib.use('Qt4Agg')
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import numpy
import os
import pyFAI
from pyFAI.gui_utils import pylab, QtGui, QtCore, uic
import sys


SIGNAL = pyFAI.gui_utils.QtCore.SIGNAL

logger = logging.getLogger("calibration_view")

#for debug
cw = None

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
        self.axes = None
        # ar is for artists: plot, images or labels...
        self.ar_data = self.ar_mask = self.ar_massif = self.ar_contour = self.ar_solidangle = self.ar_points = None
        self.data = self.massif = self.solid_angle = self.mask = None
        self.display_checks = {}
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

        self.axes = self.fig.add_subplot(111)
        self.axes.set_visible(False)
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
        self.display_checks = {"data": self.display_widget.show_data,
                               "mask": self.display_widget.show_mask,
                               "massif": self.display_widget.show_massif,
                               "points": self.display_widget.show_points,
                               "contour": self.display_widget.show_contour,
                               "solidangle": self.display_widget.show_solidangle}
        for v in self.display_checks.itervalues():
            self.connect(v, SIGNAL("stateChanged(int)"), self.toggle_show)

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

    def set_data(self, data, display=True, target="data"):
        """
        Display an array in the  data layer  
        @param data: the numpy array with the data in it
        @param display: shall the data be displayed  
        """
        self.__setattr__(target, data)
        artist = self.__getattribute__("ar_%s" % target)
        if self.axes is None:
            return
        if artist is None:
            artist = self.axes.imshow(data, zorder=self.ZORDER[target])
            self.__setattr__("ar_%s" % target, artist)
        else:
            artist.set_data(data)
        artist.set_visible(display)
        self.canvas.draw()
        if display:
            self.display_checks[target].setChecked(True)

    def calc_RGBA(self, data, target="data"):
        """
        Apply the colormap depending on the target
        
        @param data: array of floats
        @return: y,x,4 array of uint8
        
        TODO
        """
        if target == "data":
            datamin = data.min()
            showData = numpy.log1p(data - self.data.min())
#            self.ax.set_title('Log colour scale (skipping lowest/highest per mille)')
#        else:
#            showData = self.data
#            self.ax.set_title('Linear colour scale (skipping lowest/highest per mille)')

        # skip lowest and highest per mille of image values via vmin/vmax
        showMin = percentile(showData, .1)
        showMax = percentile(showData, 99.9)
        im = self.ax.imshow(showData, vmin=showMin, vmax=showMax, origin="lower", interpolation="nearest")

    def any_display(self):
        if self.axes is not None:
            display = False
            for v in self.display_checks.values():
                display = display or v.isChecked()
            self.axes.set_visible(display)
            self.canvas.draw()

    def toggle_show(self):
        if self.axes is not None:
            display = False
            for k, v in self.display_checks.iteritems():
                artist = self.__getattribute__("ar_" + k)
                if artist is not None:
                    display_artsist = v.isChecked()
                    artist.set_visible(display_artsist)
                    display = display or display_artsist
            self.axes.set_visible(display)
            self.canvas.draw()


if __name__ == "__main__":
    filename = os.path.join(os.environ["HOME"], "workspace", "pyFAI", "test",
                            "testimages", "Pilatus1M.edf")
    app = QtGui.QApplication(sys.argv)
    cw = CalibrationWindow()

    cw.show()
    cw.set_data(fabio.open(filename).data)
    det = pyFAI.detectors.Pilatus1M()
    cw.set_data(det.mask, target="mask")
    sys.exit(app.exec_())
