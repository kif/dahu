#!/usr/bin/python
import sys, os, time, json

if "TANGO_HOST" not in os.environ:
    raise RuntimeError("No TANGO_HOST defined")
import PyTango

dahu = PyTango.DeviceProxy("DAU/dahu/1")

data = {
    "output_hdf5_filename": "/nobackup/lid02gpu12/test.h5",
    "output_hdf5_dataset":"entry_000/data", 
    #"distortion_spline": "/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline", 
    "input_hdf5_filename": "/nobackup/lid02gpu11/FRELON/test_laurent_saxs_0000.h5", 
#     "output_hdf5_filename": "/data/opid02/inhouse/test_dahu/dahu.h5", 
    "input_hdf5_dataset": "/entry_0000/measurement/detector/data"
}
kwargs = "{\"x\":5}" 

for i in range(1):
    for plugin in [ "example.cube","example.square","example.sleep"]:
        pid = dahu.startJob([plugin, kwargs])
print("%s id: %i"%(plugin, pid))
print("Input: %s"% dahu.getJobInput(pid))
print("Output: %s"% dahu.getJobOutput(pid))
print("state: %s"% dahu.getJobState(pid))

def startjob(spline=None,method=None,device=None,wg=None):
    if method:
        data["method"] = method
    if wg:
        data["workgroup"] = wg
    if device:
        data["device"] = device
    if spline:
        data["spline"] = spline
    if os.path.exists(data["output_hdf5_dataset"]):
        os.unlink(data["output_hdf5_dataset"])
    print(data)
    pid = dahu.startJob(["id02.distortion", json.dumps(data)])
    print("%s id: %i"%(plugin, pid))
    print("Input: %s"% dahu.getJobInput(pid))
    print("Output: %s"% dahu.getJobOutput(pid))
    print("state: %s"% dahu.getJobState(pid))
    while dahu.getJobState(pid) not in ["success", "failure"]:
        time.sleep(1)
    print("output:")
    print(dahu.getJobOutput(pid))
    
for i in ((None,None,None,None),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline",None,None,None),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","lut",None,None),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","csr",None,None),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","lut","gpu",8),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","lut","cpu",8),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","csr","gpu",8),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","csr","cpu",8),
          ):
    startjob(*i)

dahu.collectStatistics()