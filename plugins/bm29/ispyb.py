#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

Everything to send data to Ispyb
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "03/12/2024"
__status__ = "development"
version = "0.2.2"

import logging
logger = logging.getLogger("bm29.ispyb")
import os
import shutil
import json
import tempfile
import numpy
try:
    from suds.client import Client
    from suds.transport.https import HttpAuthenticated
except ImportError:
    print("ISPyB connection will no work")
    Client = HttpAuthenticated = None
try:
    from pyicat_plus.client.main import IcatClient
except ImportError:
    print("iCat connection will no work")
    IcatClient = None

import matplotlib.pyplot
matplotlib.use("Agg")
from freesas.containers import RG_RESULT, RT_RESULT, StatsResult
from freesas.plot import kratky_plot, guinier_plot, scatter_plot, density_plot, hplc_plot


class NumpyEncoder(json.JSONEncoder):

    def default(self, obj):
        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def str_list(lst):
    "Helper function to convert list of path to smth compatible with ispyb"
    return json.dumps([{"filePath": i} for i in lst])


class IspybConnector:
    "This class is a conector to the web-service"

    def __init__(self, url, login=None, passwd=None, gallery=None, pyarch=None,
                 experiment_id=-1, run_number=-1, **kwargs):
        """Constructor of the ISPyB connections

        :param server: url of the service
        :param login: name used to authenticate
        :param passwd: password to use
        :param pyarch: folder for archiving data
        :param experiment_id: identifies the experiment
        :param run_number: identify the run in an experiment (sample/buffer localization)
        """
        self._url = url
        # try:
        self.authentication = HttpAuthenticated(username=login, password=passwd)
        self.client = Client(url, transport=self.authentication, cache=None)
        # except Exception as err:
        #     logger.error(f"Connection to the web-service {type(err)}:{err}")

        if gallery:
            self.gallery = os.path.abspath(gallery)
            if not os.path.isdir(self.gallery):
                try:
                   os.makedirs(self.gallery)
                except Exception as err:
                    logger.warning(f"Unable to create dir {self.gallery}. {type(err)}: {err}")
        else:
            logger.error("No `gallery` destination provided ... things will go wrong")
            self.gallery = tempfile.gettempdir()
        if pyarch:
            self.pyarch = os.path.abspath(pyarch)
        else:
            logger.error("No `pyarch` destination provided ... things will go wrong")

        self.experiment_id = experiment_id
        self.run_number = run_number
        

    def __repr__(self):
        return f"Ispyb connector to {self._url}"

    def send_icat(self, proposal=None, beamline=None, sample=None, dataset=None, path=None, raw=None,  data=None):
        """
        :param proposal: mx1324
        :param beamline: name of the beamline
        :param sample: sample name as registered in icat
        :param dataset: name given by BLISS
        :param path: directory name where processed data are staying
        :param raw: directory name of the raw data (not the processed ones)
        :param data: dict with all data sent to ISpyB
        """
        tmp = self.gallery.strip("/").split("/")
        idx_process = [i for i,j in enumerate(tmp) if j.lower().startswith("process")][-1]
        assert idx_process>5
        if proposal is None:
            proposal = tmp[idx_process-5]
        if beamline is None:
            beamline = tmp[idx_process-4]
        if sample is None:
            sample = tmp[idx_process-2]
        if dataset is None:
            dataset = tmp[idx_process+1]
        if path is None:
            path = os.path.dirname(self.gallery)
        if raw is None:            
            raw = os.path.abspath(self.gallery[:self.gallery.lower().index("process")])            
        
        metadata = {"definition": "SAXS",
                    "Sample_name": sample}
        for k,v in data.items():
            if isinstance(k, str) and k.startswith("SAXS_"):
                metadata[k] = v
        sample = data.get("sample")
        if sample:
            metadata["SAXS_concentration"] = str(sample.concentration)
            metadata["SAXS_code"] = sample.name
            metadata["SAXS_comments"] = sample.description
            metadata["SAXS_storage_temperature"] = str(sample.temperature_env)
            metadata["SAXS_exposure_temperature"] = str(sample.temperature)
            if sample.hplc:
                metadata["SAXS_column_type"] = sample.hplc
            #"buffer": "description of buffer, pH, ...",

        guinier = data.get("guinier")
        if guinier:
            metadata["SAXS_guinier_rg"] = f"{guinier.Rg:.1f}±{guinier.sigma_Rg:.1f}"
            metadata["SAXS_guinier_points"] = f"{guinier.start_point}-{guinier.end_point}"
            metadata["SAXS_guinier_i0"] = f"{guinier.I0:.1f}±{guinier.sigma_I0:.1f}"

        bift = data.get("bift")
        if bift:
            metadata["SAXS_rg"] =  f"{bift.Rg_avg:.1f}±{bift.Rg_std:.1f}"
            metadata["SAXS_d_max"] = f"{bift.Dmax_avg:.1f}±{bift.Dmax_std:.1f}"

        tomerge = data.get("merged")
        if tomerge:
            metadata["SAXS_frames_averaged"] = f"{tomerge[0]}-{tomerge[1]}"
        
        volume = data.get("volume")
        if volume:
            metadata["SAXS_porod_volume"] = str(volume) 
        #Other metadata one may collect ...
        metadata["SAXS_experiment_type"]= data.get("experiment_type", "UNKNOWN")
        metadata["datasetName"] = dataset
        icat_client = IcatClient(metadata_urls=["bcu-mq-01.esrf.fr:61613", "bcu-mq-02.esrf.fr:61613"])
        kwargs = {"beamline":beamline, 
                  "proposal":proposal, 
                  "dataset":dataset, 
                  "path":path, 
                  "metadata":metadata, 
                  "raw":[raw]}
        icat_client.store_processed_data(**kwargs)
        return kwargs

    def send_averaged(self, data):
        """Send this to ISPyB and backup to PyArch

        :param: data: dict to be saved in pyarch with keys:
                integer: frame index to be saved
                "avg": the averaged frame
                "merged": list of index merged
                0,1,2,3 the different indexes for individual frames.
        """
        data = data.copy()
        basename = data.pop("basename")
        discarded = []
        frames = []
        merged = data.pop("merged")

        aver_data = data.pop("avg")
        averaged = self.save_curve("avg", aver_data, basename, gallery=True)
        list_merged = list(range(*merged))
        for k, v in data.items():
            if isinstance(k, int):
                fn = self.save_curve(k, v, basename, gallery=False)
                if k in list_merged:
                    frames.append(fn)
                else:
                    discarded.append(fn)
        self.client.service.addAveraged(str(self.experiment_id),
                                        str(self.run_number),
                                        str_list(frames),
                                        str_list(discarded),
                                        str(averaged))
        sasm = numpy.vstack((aver_data.radial, aver_data.intensity, aver_data.sigma)).T
        self.scatter_plot(sasm, basename=basename)

    def _mk_filename(self, index, path, basename="frame", ext=".dat"):
        dest = os.path.join(self.pyarch, path)
        if not os.path.isdir(dest):
            try:
                os.makedirs(dest)
            except Exception as err:
                logger.error("Unable to create directory %s: %s: %s", dest, type(err), err)
        os.stat(dest) #this is to enforce the mounting of the directory
        if isinstance(index, int):
            filename = os.path.join(dest, "%s_%04d%s" % (basename, index, ext))
        else:
            filename = os.path.join(dest, "%s_%s%s" % (basename, index, ext))
        return filename

    def save_curve(self, index, integrate_result, basename="frame", gallery=False):
        """Save a  1D curve into the pyarch. Not those file do not exist outside pyarch

        :param: index: prefix or index value for
        :param: integrate_result: an IntegrationResult to be saved
        :param gallery: make a copy to the gallery
        :return: the full path of the file in pyarch
        """
        pyarch_fn = self._mk_filename(index, "1d", basename)
        sasl = numpy.vstack((integrate_result.radial, integrate_result.intensity, integrate_result.sigma))
        numpy.savetxt(pyarch_fn, sasl.T)
        if gallery:
            gallery_fn = os.path.join(self.gallery, os.path.basename(pyarch_fn))
            shutil.copyfile(pyarch_fn, gallery_fn)
        return pyarch_fn

    def save_bift(self, bift, basename="frame"):
        """Save a  IFT curve into the pyarch. Not those file do not exist outside pyarch

        :param: index: prefix or index value for
        :param: bift: an StatResults object to be saved (freesas >= 0.8.4)
        :return: the full path of the file in pyarch
        """
        pyarch_fn = self._mk_filename("BIFT", "1d", basename, ext=".out")
        gallery_fn = os.path.join(self.gallery, os.path.basename(pyarch_fn))
        bift.save(gallery_fn)
        shutil.copyfile(gallery_fn, pyarch_fn)
        return pyarch_fn

    def kratky_plot(self, sasm, guinier, basename="frame"):
        pyarch_fn = self._mk_filename("Kratky", "plot", basename, ext=".png")
        gallery_fn = os.path.join(self.gallery, os.path.basename(pyarch_fn))
        fig = kratky_plot(sasm, guinier,
                           filename=gallery_fn, img_format="png", unit="nm",
                           title="Dimensionless Kratky plot",
                           ax=None, labelsize=None, fontsize=None)
        matplotlib.pyplot.close(fig)
        shutil.copyfile(gallery_fn, pyarch_fn)
        return pyarch_fn

    def guinier_plot(self, sasm, guinier, basename="frame"):
        pyarch_fn = self._mk_filename("Guinier", "plot", basename, ext=".png")
        gallery_fn = os.path.join(self.gallery, os.path.basename(pyarch_fn))
        fig = guinier_plot(sasm, guinier, filename=gallery_fn,
                            img_format="png", unit="nm",
                            ax=None, labelsize=None, fontsize=None)
        matplotlib.pyplot.close(fig)
        shutil.copyfile(gallery_fn, pyarch_fn)
        return pyarch_fn

    def scatter_plot(self, sasm, guinier=None, ift=None, basename="frame"):
        pyarch_fn = self._mk_filename("Scattering", "plot", basename, ext=".png")
        gallery_fn = os.path.join(self.gallery, os.path.basename(pyarch_fn))
        fig = scatter_plot(sasm, guinier, ift,
                           filename=gallery_fn, img_format="png", unit="nm",
                           title="Scattering curve ",
                           ax=None, labelsize=None, fontsize=None)
        matplotlib.pyplot.close(fig)
        shutil.copyfile(gallery_fn, pyarch_fn)
        return pyarch_fn

    def density_plot(self, ift, basename="frame"):
        pyarch_fn = self._mk_filename("Density", "plot", basename, ext=".png")
        gallery_fn = os.path.join(self.gallery, os.path.basename(pyarch_fn))
        fig = density_plot(ift, filename=gallery_fn, img_format="png", unit="nm",
                     ax=None, labelsize=None, fontsize=None)
        matplotlib.pyplot.close(fig)
        shutil.copyfile(gallery_fn, pyarch_fn)
        return pyarch_fn

    def send_subtracted(self, data):
        """send the result of the subtraction to Ispyb

        :param data: a dict with all information to be saved in Ispyb
        """
        run_number = list(self.run_number)
        guinier = data.get("guinier")
        gnom = data.get("bift")
        subtracted = data.get("subtracted")
        basename = data.get("basename", "frame")
        sub = self.save_curve("subtracted", subtracted, basename, gallery=True)
        buf = self.save_curve("buffer_avg", data.get("buffer"), basename)
        individual_buffers = []
        for i, bufi in enumerate(data.get("buffers", [])):
            individual_buffers.append(self.save_curve("buffer_%d" % i, bufi, basename))

        sample = self.save_curve("sample", data.get("sample"), basename)
        if gnom is not None:
            gnomFile = self.save_bift(gnom, basename)
        sasm = numpy.vstack((subtracted.radial, subtracted.intensity, subtracted.sigma)).T
        if guinier:
            kratkyPlot = self.kratky_plot(sasm, guinier, basename)
            guinierPlot = self.guinier_plot(sasm, guinier, basename)
        else:
            kratkyPlot = ""
            guinierPlot = ""

        scatterPlot = self.scatter_plot(sasm, guinier, gnom, basename)
        if gnom is not None:
            densityPlot = self.density_plot(gnom, basename)
        else:
            densityPlot = ""

        self.client.service.addSubtraction(str(self.experiment_id),
                                           str(run_number),
                                           str(guinier.Rg if guinier else -1),
                                           str(guinier.sigma_Rg if guinier else -1),
                                           str(guinier.I0 if guinier else -1),
                                           str(guinier.sigma_I0 if guinier else -1),
                                           str(guinier.start_point if guinier else -1),
                                           str(guinier.end_point if guinier else -1),
                                           str(guinier.quality if guinier else -1),
                                           str(guinier.aggregated if guinier else -1),
                                           str(gnom.Rg_avg if gnom else -1),
                                           str(gnom.Dmax_avg if gnom else -1),
                                           str(gnom.evidence_avg if gnom else -1),
                                           str(data.get("volume", -1)),
                                           "[{'filePath': '%s'}]" % sample,  # sampleOneDimensionalFiles
                                           str_list(individual_buffers),  # bufferOneDimensionalFiles
                                           sample,  # sampleAverageFilePath,
                                           buf,  # bufferAverageFilePath,
                                           sub,  # subtractedFilePath,
                                           scatterPlot,
                                           densityPlot,
                                           guinierPlot,
                                           kratkyPlot,
                                           gnomFile if gnom else "")

    def send_hplc(self, data):
        """send the result of the subtraction to Ispyb

        :param data: a dict with all information to be saved in Ispyb
        """
        sample = data.get("sample_name", "sample")
        #gallery
        gallery = os.path.join(self.gallery, 'chromatogram.png')
        chromatogram = data.get("sum_I")
        
        if chromatogram is not None:
            fractions = data.get("merge_frames")
            if fractions is not None:
                fractions.sort()
            hplc_plot(chromatogram, fractions,
                      title=f"Chromatogram of {sample}", 
                      filename=gallery, 
                      img_format="png", )

        hdf5_file = data.get("hdf5_filename")
        filename = self._mk_filename("hplc", ".", sample, ext=".h5")
        filename = os.path.abspath(filename)
        json_file = os.path.splitext(filename)[0] + ".json"
        with open(json_file, mode="w") as w:
            w.write(json.dumps(data, indent=2, cls=NumpyEncoder))
        shutil.copyfile(hdf5_file, filename)
        self.client.service.storeHPLC(str(self.experiment_id),
                                      filename,
                                      json_file)

