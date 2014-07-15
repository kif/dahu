#!/usr/bin/python
import sys, os, time, json

if "TANGO_HOST" not in os.environ:
    raise RuntimeError("No TANGO_HOST defined")
import PyTango

dahu = PyTango.DeviceProxy("DAU/dahu/1")



kwargs = "{\"x\":5}"

for i in range(1):
    for plugin in [ "example.cube", "example.square", "example.sleep"]:
        pid = dahu.startJob([plugin, kwargs])
print("%s id: %i" % (plugin, pid))
print("Input: %s" % dahu.getJobInput(pid))
print("Output: %s" % dahu.getJobOutput(pid))
print("state: %s" % dahu.getJobState(pid))
# dahu.collectStatistics()
# print(dahu.getStatistics())


def startjob(method, device, wg):
    data_input = json.loads(open("demo_distortion.json").read())
    data_input["method"] = method
    data_input["device"] = device
    data_input["workgroup"] = wg
    data = json.dumps(data_input)
    print(data)
    pid = dahu.startJob(["id02.distortion", data])
    print("%s id: %i" % (plugin, pid))
    print("Input: %s" % dahu.getJobInput(pid))
    print("Output: %s" % dahu.getJobOutput(pid))
    print("state: %s" % dahu.getJobState(pid))
    while dahu.getJobState(pid) not in ["success", "failure"]:
        time.sleep(1)

# for device in (None,"CPU", "GPU"):
#     for method in ("csr","lut"):
#         for wg in [1,2,4,8,16,32]:
#             if (method != "CSR" or device == "None")and wg > 1:
#                 continue
#             else:
#                 startjob(method,device,wg)

data = {
        "hdf5_filename":"/nobackup/lid02gpu11/metadata/test.h5",
        "entry": "entry",
        "instrument":"id02",
        "c216":"id02/c216/0",
        "HS32F": [1e-06, 1, 7763480, 8176290, 342239, 967341, 5541980, 1739160, 2753.61, 1351920, 140000000, 16719500, 1, 0.000995868, 0.000995868, 1],
        "HS32Z": [0, 0, 383.55, 126.4, 6582.1, 6973.6, 162.95, 0, 221.2, 207.8, 315.26, 145.11, 323.76, 170, 170, 228.4],
        "HS32N": ["TIME", "AUX1", "PIN41", "PIN42", "PIN5", "PIN6", "PIN7", "PIN8", "PIN1", "PIN2", "PIN3", "PIN4", "AUX2", "THC1", "THC2", "PRESS"],
        "HSI0": 12,
        "HSI1": 7,
        "HSTime": 1,
        "HMStartEpoch": 1405087717.12159,
        "HMStartTime": "2014-07-11T16:08:37.121591+0200",
        "info": {"DetectorInfo":"VxH:detbox=14952.235x0.000x1.000,dettab=-62.000x-245.000",
                 "ExperimentInfo":"0",
                 "MachineInfo": "Ie=183.27mA,u35u=100.000mm/0.001mm,u21m=100.000mm/0.000mm,u21d=100.000mm/-0.000mm",
                 "MirrorInfo": "rz=-3.600mrad,ty=0.455mm,ty1=2.075mm,ty2=-1.165mm,tz1=-0.030mm,tz2=-0.090mm,mir2rz=2.000mrad",
                 "OpticsInfo": "egy=12460.0eV,theta=9.132deg,hgt=11.7mm,tilt=4.440deg,tra=1.963mm",
                 "ProposalInfo": 0,
                 "StationInfo": "ID02"
                 }
        }

pid = dahu.startJob(["id02.metadata", json.dumps(data)])
print("%s id: %i" % (plugin, pid))
print("Input: %s" % dahu.getJobInput(pid))
print("Output: %s" % dahu.getJobOutput(pid))
print("state: %s" % dahu.getJobState(pid))
while dahu.getJobState(pid) not in ["success", "failure"]:
    time.sleep(1)

# dahu.collectStatistics()
# time.sleep(5)
# print(dahu.getStatistics())
print("Input: %s" % dahu.getJobInput(pid))
print("Output: %s" % dahu.getJobOutput(pid))
print("state: %s" % dahu.getJobState(pid))
print("error: %s" % dahu.getJobError(pid))
