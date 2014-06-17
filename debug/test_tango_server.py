#!/usr/bin/python
import sys, os, time

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



data_input = open("demo_distortion.json").read()
pid = dahu.startJob(["id02.distortion", data_input])
print("%s id: %i"%(plugin, pid))
print("Input: %s"% dahu.getJobInput(pid))
print("Output: %s"% dahu.getJobOutput(pid))
print("state: %s"% dahu.getJobState(pid))
while dahu.getJobState(pid) not in ["success", "failure"]:
    time.sleep(1)

dahu.collectStatistics()
time.sleep(1)
print(dahu.getStatistics())
print("Input: %s"% dahu.getJobInput(pid))
print("Output: %s"% dahu.getJobOutput(pid))
print("state: %s"% dahu.getJobState(pid))
