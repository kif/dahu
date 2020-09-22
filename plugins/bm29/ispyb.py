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
import numpy    
from suds.client            import Client
from suds.transport.https   import HttpAuthenticated
from freesas.collections import RG_RESULT

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
            
        for k,v in data.items():
            if isinstance(k, int):
                fn = self.save_curve(k, v, basename)
                if k in merged:
                    frames.append(fn)
                else:
                    discarded.append(fn)
        self.client.service.addAveraged(str(self.measurement_id),
                                        str(self.collection_id),
                                        str(frames),
                                        str(discarded),
                                        str(averaged))
    
    def save_curve(self, index, integrate_result, basename="frame"):
        """Save a  1D curve into the pyarch. Not those file do not exist outside pyarch
        
        :param: index: prefix or index value for 
        :param: integrate_result: an IntegrationResult to be saved.  
        :return: the full path of the file in pyarch 
        """
        dest = os.path.join(self.pyarch, "1d")
        if not os.path.isdir(dest):
            os.makedirs(dest)
        if isinstance(index, int):
            filename = os.path.join(dest, "%s_%04d.dat"%(basename,index))
        else:
            filename = os.path.join(dest, "%s_%s.dat"%(basename, index))
        sasl = numpy.vstack((integrate_result.radial, integrate_result.intensity, integrate_result.sigma))
        numpy.savetxt(filename, sasl.T)
        return filename
        
    def send_subtracted(self, data):
        """send the result of the subtraction to Ispyb
        
        :param data: a dict with all information to be saved in Ispyb
        """ 
        guinier = data.get("guinier", *([-1.]*8))
        gnom = data.get("bift", *([-1.]*8))
        subtracted = data.get("subtracted")
        sub = self.save_curve("subtracted", subtracted, data.get("basename", "frame"))
        sasm = numpy.vstack((subtracted.radius, subtracted.intensity, subtracted.sigma))
        
        self.client.service.addSubtraction(str(self.measurement_id),
                                           str(guinier.Rg),
                                           str(guinier.sigma_Rg),
                                           str(guinier.I0),
                                           str(guinier.sigma_I0),
                                           str(guinier.start_point),
                                           str(guinier.end_point),
                                           str(guinier.quality),
                                           str(guinier.aggregated),
                                           str(gnom.Rg_avg),
                                           str(gnom.Dmax_avg),
                                           str(gnom.logP_avg),
                                           str(data["volume"]),
                                           str(sampleAvgOneDimensionalFiles),
                                           str(bufferAvgOneDimensionalFiles),
                                           self.averageSample,                     #sampleAverageFilePath,
                                           self.bestBuffer,                        #bufferAverageFilePath,
                                           sub,                                    #subtractedFilePath,
                                           self.scatterPlot,                       #experimentalDataPlotFilePath,
                                           self.densityPlot,                       #densityPlotFilePath,
                                           self.guinierPlot,                       #guinierPlotFilePath,
                                           self.kratkyPlot,                        #kratkyPlotFilePath,
                                           self.gnomFile)                          #gnomOutputFilePath
