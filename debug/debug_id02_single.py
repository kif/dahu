#!/usr/bin/python
import sys, os, time, json

if "TANGO_HOST" not in os.environ:
    raise RuntimeError("No TANGO_HOST defined")
import PyTango

dahu = PyTango.DeviceProxy("DAU/dahu/1")
plugin = "id02.SingleDetector"
data = {"DetectorName": 'Rayonix',
        "image_file": "/nobackup/lid02gpu11/FRELON/test_laurent_saxs_0000.h5",
        #"entry": "entry_0000"
        #"hdf5": "/entry_0000/id02/data
        "output_dir": "/nobackup/lid02gpu12",
        "PSize_1": 2.4e-05,
        "PSize_2": 2.4e-05,
        "BSize_1":1,
        "BSize_2":1,
        "Center_1": 512,
        "Center_2": 512,
        "DDummy": 1,
        "SampleDistance":14.9522,
        "c216_filename": "/nobackup/lid02gpu11/metadata/test.h5",
        "WaveLength": 9.95058e-11,
        "Dummy":-10,
        "output_dir": "/nobackup/lid02gpu12/output",
        #              "do_dark":false,
        #              "do_flat":false,
        "npt2_azim": 360,
        "npt2_rad" : 500,
        "npt1_rad" : 1000,
        "to_save": ["raw", "azim", "ave"],
        }
print dahu.initPlugin(plugin)
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
