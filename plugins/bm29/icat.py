#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

Everything to send data to iCat, the data catalogue
 
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "21/02/2025" 
__status__ = "development"
version = "0.3.0"


import os
import tempfile
import logging
logger = logging.getLogger(__name__)
try:
    from pyicat_plus.client.main import IcatClient
except ImportError:
    logger.error("iCat connection will no work")
    IcatClient = None

def _ensure_gallery(gallery):
    if gallery:
        gallery = os.path.abspath(gallery)
        if not os.path.isdir(gallery):
            try:
                os.makedirs(gallery)
            except Exception as err:
                logger.warning(f"Unable to create dir {gallery}. {type(err)}: {err}")
    else:
        logger.error("No `gallery` destination provided ... things will go wrong")
        gallery = tempfile.gettempdir()
    return gallery


def send_icat(proposal=None, beamline=None, sample=None, dataset=None, path=None, raw=None,  data=None, gallery=None, metadata=None):
    """Send some data to icat, the data-catalogue
    
    :param proposal: mx1324
    :param beamline: name of the beamline
    :param sample: sample name as registered in icat
    :param dataset: name given by BLISS
    :param path: directory name where processed data are staying
    :param raw: list of directory name of the raw data (not the processed ones)
    :param data: dict with all data sent to iCat
    :param gallery: path with the gallery directory
    :param metadata: dict with additionnal metadata (could be overwritten by this function
    :return: data sent to icat as a dict
    """
    gallery = _ensure_gallery(gallery)
    tmp = gallery.strip("/").split("/")
    idx_process = [i for i,j in enumerate(tmp) if j.lower().startswith("process")][-1]
    if tmp[idx_process] == "processed":
        assert idx_process>=6
        if proposal is None:
            proposal = tmp[idx_process-6]
        if beamline is None:
            beamline = tmp[idx_process-5]
        if sample is None:
            sample = tmp[idx_process-2]
        if dataset is None:
            dataset = tmp[idx_process+1]
        if path is None:
            path = os.path.dirname(gallery)
        if raw is None:            
            raw = os.path.abspath(gallery[:gallery.lower().index("process")])
    elif tmp[idx_process] == "PROCESSED_DATA":           
        if proposal is None:
            proposal = tmp[idx_process-3]
        if beamline is None:
            beamline = tmp[idx_process-2]
        if sample is None:
            sample = tmp[idx_process+1]
        if dataset is None:
            dataset = tmp[idx_process+2]
        if path is None:
            path = os.path.dirname(gallery)
        if raw is None:            
            raw = os.path.dirname(os.path.dirname(os.path.abspath(gallery.replace("PROCESSED_DATA", "RAW_DATA"))))
    else:
        logger.error("Unrecognized path layout")
    
    if metadata is None:
        metadata = {}
    metadata["definition"] = "SAXS",
   # metadata["Sample_name"] = sample
    
    for k,v in data.items():
        if isinstance(k, str) and k.startswith("SAXS_"):
            metadata[k] = v
    sample = data.get("sample", sample)
    if sample is not None:
        if  isinstance(sample, str):
            metadata["Sample_name"] = metadata["SAXS_code"] = sample
        else:
            metadata["SAXS_concentration"] = str(sample.concentration)
            metadata["Sample_name"] = metadata["SAXS_code"] = sample.name
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
        metadata["SAXS_chi2r"] = str(bift.chi2r_avg)
        metadata["SAXS_chi2r_error"] = str(bift.chi2r_std)

    tomerge = data.get("merged")
    if tomerge:
        metadata["SAXS_frames_averaged"] = f"{tomerge[0]}-{tomerge[1]}"
    
    volume = data.get("volume")
    if volume:
        metadata["SAXS_porod_volume"] = str(volume) 
    rti = data.get("rti")
    if rti:
        "Vc sigma_Vc Qr sigma_Qr mass sigma_mass"
        metadata["SAXS_vc"] = str(rti.Vc)
        metadata["SAXS_vc_error"] = str(rti.sigma_Vc)
        metadata["SAXS_mass"] = str(rti.mass)
        metadata["SAXS_mass_error"] = str(rti.sigma_mass)
        
    if not isinstance(raw, list):
        raw = [raw]
    #Other metadata one may collect ...
    metadata["SAXS_experiment_type"]= data.get("experiment_type", "UNKNOWN")
    metadata["datasetName"] = dataset
    icat_client = IcatClient(metadata_urls=["bcu-mq-01.esrf.fr:61613", "bcu-mq-02.esrf.fr:61613"])
    kwargs = {"beamline":beamline, 
              "proposal":proposal, 
              "dataset":dataset, 
              "path":path, 
              "metadata":metadata, 
              "raw":raw}
    #print(kwargs)
    icat_client.store_processed_data(**kwargs)
    return kwargs
