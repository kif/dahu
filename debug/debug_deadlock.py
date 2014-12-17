#!/usr/bin/python

import os
import PyTango, json

saxs = {
    "npt2_rad": 1000, 
    "dark_filter_quantil_lower": 0.1, 
    "dark_filter_quantil_upper": 0.9, 
    "regrouping_mask_filename": "mask-05m.edf", 
    "npt2_azim": 360, 
    "job_id": 2, 
    "dark_filename": "/nobackup/lid02gpu12/nl16b_saxs_02229_dark.h5", 
    "to_save": "raw sub flat norm azim ave", 
    "dark_filter": "quantil", 
    "output_dir": "/mntdirect/_data_opid02_inhouse/com/20141216d/cor", 
    "image_file": "/nobackup/lid02gpu12/nl16b_saxs_02229.h5", 
    "npt1_rad": 1000, 
    "dark_filter_quantil": 0.5, 
    "DetectorName": "saxs", 
    "plugin_name": "id02.singledetector", 
    "flat_filename": "/data/opid02/archive/setup/spatcorr-files/saxs/flat_saxs_2x2.edf", 
    "distortion_filename": "/data/opid02/archive/setup/spatcorr-files/saxs/SpatCorrRayonix_2b2.dat", 
    "c216_filename": "/nobackup/lid02gpu12/nl16b_scalers_02229.h5"}

waxs = {
    "npt2_rad": 1000, 
    "dark_filter_quantil_lower": 0.1, 
    "dark_filter_quantil_upper": 0.9, 
    "regrouping_mask_filename": "mask-waxs.edf", 
    "npt2_azim": 360, 
    "job_id": 3, 
    "dark_filename": "/nobackup/lid02gpu11/nl16b_waxs_02229_dark.h5", 
    "to_save": "raw sub flat norm azim ave", 
    "dark_filter": "quantil", 
    "output_dir": "/mntdirect/_data_opid02_inhouse/com/20141216d/cor", 
    "image_file": "/nobackup/lid02gpu11/nl16b_waxs_02229.h5", 
    "npt1_rad": 1000, 
    "dark_filter_quantil": 0.5, 
    "DetectorName": "waxs", 
    "plugin_name": "id02.singledetector", 
    "c216_filename": "/nobackup/lid02gpu12/nl16b_scalers_02229.h5"
}

dahu = PyTango.DeviceProxy("DAU/dahu/1")
js = dahu.startJob(["id02.singledetector",json.dumps(saxs)])
jw = dahu.startJob(["id02.singledetector",json.dumps(waxs)])

dahu.collectStatistics()
print(dahu.statisticsCollected)