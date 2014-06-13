#!/usr/bin/python
import sys, os, time

if "TANGO_HOST" not in os.environ:
    raise RuntimeError("No TANGO_HOST defined")
import PyTango

dahu = PyTango.DeviceProxy("DAU/dahu/1")



kwargs = "{\"x\":5}" 

for i in range(1000):
    for plugin in [ "example.cube","example.square","example.sleep"]:
        pid = dahu.startJob([plugin, kwargs])
        print("%s id: %i"%(plugin, pid))
        print("Input: %s"% dahu.getJobInput(pid))
        print("Output: %s"% dahu.getJobOutput(pid))
        print("state: %s"% dahu.getJobState(pid))