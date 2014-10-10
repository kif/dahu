#!/usr/bin/python
import sys, os, time, json, glob

if "TANGO_HOST" not in os.environ:
    raise RuntimeError("No TANGO_HOST defined")
import PyTango

dahu = PyTango.DeviceProxy("DAU/dahu/1")

def status(pid):
    try:
        res = dahu.getJobState(pid)
    except:
        print("Communication failure with Dahu via Tango")
    else:
        return res


data = {
    #"output_hdf5_filename_base": "/nobackup/lid02gpu12/test",
    "output_hdf5_filename_base": "/data/opid02/inhouse/test_dahu/test",
    "output_hdf5_dataset":"entry_000/data", 
    #"distortion_spline": "/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline", 
    "input_hdf5_filename": "/nobackup/lid02gpu11/FRELON/test_laurent_saxs_0000.h5", 
 
    "input_hdf5_dataset": "/entry_0000/measurement/detector/data"
}
kwargs = "{\"x\":5}" 

for i in range(1):
    for plugin in [ "example.cube","example.square","example.sleep"]:
        pid = dahu.startJob([plugin, kwargs])
print("%s id: %i"%(plugin, pid))
print("Input: %s"% dahu.getJobInput(pid))
print("Output: %s"% dahu.getJobOutput(pid))
print("state: %s"% status(pid))

    
def startjob(spline=None,method=None,device=None,wg=None):
    if method:
        data["method"] = method
    if wg:
        data["workgroup"] = wg
    if device:
        data["device"] = device
    if spline:
        data["distortion_spline"] = spline
    data["output_hdf5_filename"] = "%s_%04i.h5"%(data["output_hdf5_filename_base"],len(glob.glob(data["output_hdf5_filename_base"]+"*")))
    print(data)
    pid = dahu.startJob(["id02.distortion", json.dumps(data)])
    print("%s id: %i"%(plugin, pid))
    print("Input: %s"% dahu.getJobInput(pid))
    print("Output: %s"% dahu.getJobOutput(pid))
    print("state: %s"% status(pid))
    while status(pid) not in ["success", "failure"]:
        time.sleep(1)
    print("output:")
    output =dahu.getJobOutput(pid) 
    print(output)
    output
    print("Average write speed: %f MB/s"%(os.path.getsize(data["output_hdf5_filename"])/json.loads(output)["job_runtime"]/2**20))
    
for i in ((None,None,None,None),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline",None,None,None),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","lut",None,None),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","csr",None,None),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","lut","gpu",8),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","lut","cpu",8),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","csr","gpu",8),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","csr","cpu",8),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","lut","gpu",32),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","lut","cpu",32),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","csr","gpu",32),
          ("/users/kieffer/workspace/pyFAI/test/testimages/frelon.spline","csr","cpu",32),

          ):
    startjob(*i)

dahu.collectStatistics()