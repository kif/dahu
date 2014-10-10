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
