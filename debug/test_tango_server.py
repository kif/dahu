#!/usr/bin/python
import sys, os, time, json

if "TANGO_HOST" not in os.environ:
    raise RuntimeError("No TANGO_HOST defined")
import PyTango

dahu = PyTango.DeviceProxy("DAU/dahu/1")



kwargs = "{\"x\":5}" 

for i in range(1):
    for plugin in [ "example.cube","example.square","example.sleep"]:
        pid = dahu.startJob([plugin, kwargs])
print("%s id: %i"%(plugin, pid))
print("Input: %s"% dahu.getJobInput(pid))
print("Output: %s"% dahu.getJobOutput(pid))
print("state: %s"% dahu.getJobState(pid))
dahu.collectStatistics()
print(dahu.getStatistics())


def startjob(method,device,wg):
    data_input = json.loads(open("demo_distortion.json").read())
    data_input["method"] = method
    data_input["device"] = device
    data_input["workgroup"] = wg
    data = json.dumps(data_input)
    print(data)
    pid = dahu.startJob(["id02.distortion", data])
    print("%s id: %i"%(plugin, pid))
    print("Input: %s"% dahu.getJobInput(pid))
    print("Output: %s"% dahu.getJobOutput(pid))
    print("state: %s"% dahu.getJobState(pid))
    while dahu.getJobState(pid) not in ["success", "failure"]:
        time.sleep(1)
    
for device in (None,"CPU", "GPU"):
    for method in ("csr","lut"):
        for wg in [1,2,4,8,16,32]:
            if (method != "CSR" or device == "None")and wg > 1:
                continue
            else:
                startjob(method,device,wg)

dahu.collectStatistics()
time.sleep(1)
print(dahu.getStatistics())
