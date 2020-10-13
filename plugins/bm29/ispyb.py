#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

Everything to send data to Ispyb
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "22229/2020"
__status__ = "development"
version = "0.0.2"

import logging
logger = logging.getLogger("bm29.ispyb")
import os
import shutil
import json
import numpy
from suds.client            import Client
from suds.transport.https   import HttpAuthenticated
from freesas.collections import RG_RESULT, RT_RESULT, StatsResult
from freesas.plot import kratky_plot, guinier_plot, scatter_plot, density_plot


def str_list(lst):
    "Helper function to convert list of path to smth compatible with ispyb"
    return json.dumps([{"filePath": i} for i in lst])


class IspybConnector:
    "This class is a conector to the web-service"

    def __init__(self, url, login=None, passwd=None, pyarch=None, collection_id=-1, measurement_id=-1):
        """Constructor of the ISPyB connections
        
        :param server: url of the service
        :param login: name used to authenticate
        :param passwd: password to use
        :param pyarch: folder for archiving data
        :param collection_id: identifier for the collection
        :param measurement_id: identifier for the measurement
        """
        self.authentication = HttpAuthenticated(username=login, password=passwd)
        self.client = Client(url, transport=self.authentication, cache=None)
        if pyarch:
            self.pyarch = os.path.abspath(pyarch)
        else:
            logger.error("No `pyarch` destination provided ... things will go wrong")

        self.collection_id = collection_id
        self.measurement_id = measurement_id

    def send_averaged(self, data):
        """Send this to ISPyB and backup to PyArch
        
        :param: data: dict to be saved in pyarch with keys:
                integer: frame index to be saved
                "avg": the averaged frame
                "merged": list of index merged
                0,1,2,3 the different indexes for individual frames.     
        """
        basename = data.pop("basename")
        discarded = []
        frames = []
        merged = data.pop("merged")

        aver_data = data.pop("avg")
        averaged = self.save_curve("avg", aver_data, basename)

        for k, v in data.items():
            if isinstance(k, int):
                fn = self.save_curve(k, v, basename)
                if k in merged:
                    frames.append(fn)
                else:
                    discarded.append(fn)
        self.client.service.addAveraged(str(self.measurement_id),
                                        str(self.collection_id),
                                        str_list(frames),
                                        str_list(discarded),
                                        str(averaged))

    def _mk_filename(self, index, path, basename="frame", ext=".dat"):
        dest = os.path.join(self.pyarch, path)
        if not os.path.isdir(dest):
            os.makedirs(dest)
        if isinstance(index, int):
            filename = os.path.join(dest, "%s_%04d%s" % (basename, index, ext))
        else:
            filename = os.path.join(dest, "%s_%s%s" % (basename, index, ext))
        return filename

    def save_curve(self, index, integrate_result, basename="frame"):
        """Save a  1D curve into the pyarch. Not those file do not exist outside pyarch
        
        :param: index: prefix or index value for 
        :param: integrate_result: an IntegrationResult to be saved.  
        :return: the full path of the file in pyarch 
        """
        filename = self._mk_filename(index, "1d", basename)
        sasl = numpy.vstack((integrate_result.radial, integrate_result.intensity, integrate_result.sigma))
        numpy.savetxt(filename, sasl.T)
        return filename

    def save_bift(self, bift, basename="frame"):
        """Save a  IFT curve into the pyarch. Not those file do not exist outside pyarch
        
        :param: index: prefix or index value for 
        :param: bift: an StatResults object to be saved (freesas >= 0.8.4).  
        :return: the full path of the file in pyarch 
        """
        filename = self._mk_filename("BIFT", "plot", basename, ext=".out")
        bift.save(filename)
        return filename

    def kratky_plot(self, sasm, guinier, basename="frame"):
        filename = self._mk_filename("Kratky", "plot", basename, ext=".png")
        kratky_plot(sasm, guinier,
                    filename=filename, format="png", unit="nm",
                    title="Dimensionless Kratky plot",
                    ax=None, labelsize=None, fontsize=None)
        return filename
    
    def guinier_plot(self, sasm, guinier, basename="frame"):
        filename = self._mk_filename("Guinier", "plot", basename, ext=".png")
        guinier_plot(sasm, guinier, filename=filename,
                 format="png", unit="nm", 
                 ax=None, labelsize=None, fontsize=None)
        return filename
    
    def scatter_plot(self, sasm, guinier, ift, basename="frame"):
        filename = self._mk_filename("Scattering", "plot", basename, ext=".png")
        scatter_plot(sasm, guinier, ift,
                 filename=filename, format="png", unit="nm",
                 title="Scattering curve ",
                 ax=None, labelsize=None, fontsize=None)
        return filename
    
    def density_plot(self, ift, basename="frame"):
        filename = self._mk_filename("Density", basename, ext=".png")
        density_plot(ift, filename=filename, format="png", unit="nm",
                     ax=None, labelsize=None, fontsize=None)
        return filename
        
    def send_subtracted(self, data):
        """send the result of the subtraction to Ispyb
        
        :param data: a dict with all information to be saved in Ispyb
        """
        guinier = data.get("guinier", RG_RESULT(*([-1.] * 8)))
        rti = data.get("rti", RT_RESULT(*([-1.] * 6)))
        volume = data.get("volume", -1) 
        gnom = data.get("bift", None)
        subtracted = data.get("subtracted")
        basename = data.get("basename", "frame")
        sub = self.save_curve("subtracted", subtracted, basename)
        buf = self.save_curve("buffer", data.get("buffer"), basename)
        sample = self.save_curve("sample", data.get("sample"), basename)
        if gnom is not None: 
            gnomFile = self.save_bift(gnom, basename)
        sasm = numpy.vstack((subtracted.radial, subtracted.intensity, subtracted.sigma)).T
        kratkyPlot = self.kratky_plot(sasm, guinier, basename)
        guinierPlot = self.guinier_plot(sasm, guinier, basename)
        scatterPlot = self.scatter_plot(sasm, guinier, gnom, basename)
        if gnom is not None:
            densityPlot = self.density_plot(gnom, basename)
        else:
            densityPlot = None

        self.client.service.addSubtraction(str(self.measurement_id),
                                           str(guinier.Rg),
                                           str(guinier.sigma_Rg),
                                           str(guinier.I0),
                                           str(guinier.sigma_I0),
                                           str(guinier.start_point),
                                           str(guinier.end_point),
                                           str(guinier.quality),
                                           str(guinier.aggregated),
                                           str(gnom.Rg_avg if gnom else -1),
                                           str(gnom.Dmax_avg if gnom else -1),
                                           str(gnom.evidence_avg if gnom else -1),
                                           str(volume),
                                           "[{'filePath': '%s'}]"%sample,  ##sampleOneDimensionalFiles
                                           "[{'filePath': '%s'}]"%buf,  ##bufferOneDimensionalFiles
                                           sample,  ##sampleAverageFilePath,
                                           buf,  ##bufferAverageFilePath,
                                           sub,                                    #subtractedFilePath,
                                           scatterPlot,
                                           densityPlot if gnom else "",
                                           guinierPlot,
                                           kratkyPlot, 
                                           gnomFile if gnom else "")          

