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
HDF5 Dialog

Simple Dialog box to display the structure of a HDF5 file

"""

__author__ = u"Jérôme Kieffer"
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "GPLv3+"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "03/09/2020"
__status__ = "development"

import logging
import os, sys
import os.path as op
from pyFAI.gui_utils import QtGui, QtCore
try:
    from PyMca5.PyMcaGUI.io.hdf5 import HDF5Widget
except ImportError:
    import HDF5Widget #use the locally provided


class HDF5Dialog(QtGui.QDialog):
    '''
    Like QFileDialog but for HDF5 data structures...
    '''

    def __init__(self, parent=None, filename=None):
        '''
        Constructor
        '''
        QtGui.QDialog.__init__(self, parent)
        self.filename = str(filename)
        self.setWindowTitle("HDF5 browser: %s" % self.filename)
        self.setModal(True)
        layout = QtGui.QVBoxLayout(self)
        self.fileModel = HDF5Widget.FileModel()
        self.fileView = HDF5Widget.HDF5Widget(self.fileModel)
        self.fileModel.openFile(self.filename)
        layout.addWidget(self.fileView)

        # OK and Cancel buttons
        self.buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
            QtCore.Qt.Horizontal, self)
        layout.addWidget(self.buttons)
        self.paths = None
        self.connect(self.buttons, QtCore.SIGNAL("accepted()"), self.accept)
        self.connect(self.buttons, QtCore.SIGNAL("rejected()"), self.close)

    def get_path(self):
        """
        @return: selected path as a string
        """
        modelIndexList = self.fileView.selectedIndexes()
        for modelIndex in modelIndexList:
            item = self.fileModel.getProxyFromIndex(modelIndex)
            return str(item.name)

    def accept(self):
        self.paths = self.get_path()
        self.close()


    @staticmethod
    def getPath(filename=None):
        """

        @return accepted, [list of entries]
        """
        dialog = HDF5Dialog(None, filename)
        dialog.exec_()
        return dialog.paths


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    f = QtGui.QFileDialog.getOpenFileName()
    print(f)
    d = HDF5Dialog.getPath(f)
    print(d)
