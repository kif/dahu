#!/usr/bin/python
import sys, os, time, json

if "TANGO_HOST" not in os.environ:
    raise RuntimeError("No TANGO_HOST defined")
import PyTango

dahu = PyTango.DeviceProxy("DAU/dahu/1")
plugin = "id02.filter"
data = { "image_file": "/nobackup/lid02gpu11/FRELON/test_laurent_saxs_0000.h5",
              #"entry": "entry_0000"
              "output_format": "edf",
              "output_file": "/nobackup/lid02gpu12/dark.edf",
              "filter": "median",
              #"cutoff": 0
              #threshold:0,
              #"format": "edf",
              #"dark_file": filename,
              #"flat_file": filename,
#              "do_dark":false,
#              "do_flat":false,
              }
print(dahu.initPlugin(plugin))
pid = dahu.startJob([plugin, json.dumps(data)])
print("%s id: %i" % (plugin, pid))
print("Input: %s" % dahu.getJobInput(pid))
print("Output: %s" % dahu.getJobOutput(pid))
print("state: %s" % dahu.getJobState(pid))
while dahu.getJobState(pid) not in ["success", "failure"]:
    time.sleep(1)
print("output:")
print(dahu.getJobOutput(pid))
if dahu.getJobState(pid) == "failure":
    print("Error:")
    print(dahu.getJobError(pid))

data = {
    "dark_filter_loq": 0.9,
    "npt2_rad": 1000,
    "dark_filter_quantil_lower": 0.1,
    "dark_filter_quantil_upper": 0.9,
    "regrouping_mask_filename": "mask-10m.edf",
    "npt2_azim": 360,
    "job_id": 66,
    "dark_filename": "/mntdirect/_data_opid02_inhouse/com/20141104/nk04_saxs_-0001_dark.h5",
    "to_save": "[\"azim\",\"ave\"]",
    "dark_filter": "quantil",
    "output_dir": "/mntdirect/_data_opid02_inhouse/com/20141104/cor",
    "image_file": "/mntdirect/_data_opid02_inhouse/com/20141104/nk04_saxs_-0001.h5",
    "npt1_rad": 1000,
    "dark_filter_quantil": 0.5,
    "DetectorName": "saxs",
    "plugin_name": "id02.singledetector",
    "flat_filename": "/data/opid02/archive/setup/spatcorr-files/saxs/flat_saxs_1x1.h5",
    "dark_filter_pq": 0.5,
    "distortion_filename": "/data/opid02/archive/setup/spatcorr-files/saxs/spline_saxs_1x1.dat",
    "c216_filename": "/mntdirect/_data_opid02_inhouse/com/20141104/nk04_scalers_-0001.h5"
}
