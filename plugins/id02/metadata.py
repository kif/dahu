"""
Metadata collection pluging which reads a C216 component using Tango 

"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "04/06/2020"
__status__ = "development"
__version__ = "0.9.1"

import os
import posixpath
import json
import numpy
from dahu.plugin import Plugin
from dahu import version as dahu_version
from dahu.utils import get_isotime, fully_qualified_name
import PyTango
import logging
logger = logging.getLogger("id02.metadata")
from .common import StringTypes, Nexus


def preproc(**d):
    """Take a dict as input and forms a metadata structure as output
    @param: any dict
    """
    dd = d.copy()
    if "job_id" in dd:
        dd.pop("job_id")
    list_f = []
    list_n = []
    list_z = []
    HS32Len = dd.get("HS32Len", 16)
    HS32Depth = dd.get("HS32Depth", 32)
    HSI0Factor = dd.get("HSI0Factor", 1)
    HSI1Factor = dd.get("HSI1Factor", 1)
    # "0.005 s"
    if "ShutterOpeningTime" in dd:
        value = dd["ShutterOpeningTime"]
        if type(value) in StringTypes:
            ShutterOpeningTime = float(value.split()[0])
        else:
            ShutterOpeningTime = float(value)
    else:
        ShutterOpeningTime = 0
    if "ShutterClosingTime" in dd:
        value = dd["ShutterClosingTime"]
        if type(value) in StringTypes:
            ShutterClosingTime = float(value.split()[0])
        else:
            ShutterClosingTime = float(value)
    else:
        ShutterClosingTime = 0
    for ind in map(lambda x: 'HS32F' + '{0:02d}'.format(x), range(1, HS32Len + 1)):
            list_f.append(float(dd[ind]))
    for ind in map(lambda x: 'HS32N' + '{0:02d}'.format(x), range(1, HS32Len + 1)):
            list_n.append(dd[ind])
    for ind in map(lambda x: 'HS32Z' + '{0:02d}'.format(x), range(1, HS32Len + 1)):
            list_z.append(float(dd[ind]))

    info_dir = {}
    for info_ind in dd:
        if info_ind[0:2].find('HS') == 0:
            continue
        elif info_ind[0:2].find('HM') == 0:
            continue
        else:
            info_dir[info_ind] = dd[info_ind]

    final_dir = {"HS32Len": HS32Len,
                 "HS32Depth": HS32Depth,
                 "HSI0Factor": HSI0Factor,
                 "HSI1Factor": HSI1Factor,
                 "ShutterOpeningTime": ShutterOpeningTime,
                 "ShutterClosingTime": ShutterClosingTime,
                 'instrument': 'id02',
                 'c216': 'id02/c216/0',
                 'HS32F': list_f,
                 'HS32Z': list_z,
                 'HS32N': list_n,
                 'Info': info_dir}
    for key in ['HMStartEpoch', 'HMStartTime', "hdf5_filename", "entry", "HSTime", "HSI0", "HSI1"]:
        if key in dd:
            final_dir[key] = dd[key]
    return final_dir


class Metadata(Plugin):
    """Plugin in charge of retrieving all metadata for ID02 and storing them into a HDF5 file

    TODO: rewrite using Nexus class from pyFAI: shorter code

    NOTA: pin number are 1-based (I0, I1, time)

    Structure of the input data:
input = {
        "hdf5_filename":"/nobackup/lid02gpu11/metadata/test.h5",
        "entry": "entry",
        "instrument":"ESRF-ID02",
        "c216":"id02/c216/0",
        "HS32F": [1e-06, 1, 7763480, 8176290, 342239, 967341, 5541980, 1739160, 2753.61, 1351920, 140000000, 16719500, 1, 0.000995868, 0.000995868, 1],
        "HS32Z": [0, 0, 383.55, 126.4, 6582.1, 6973.6, 162.95, 0, 221.2, 207.8, 315.26, 145.11, 323.76, 170, 170, 228.4],
        "HS32N": ["TIME", "AUX1", "PIN41", "PIN42", "PIN5", "PIN6", "PIN7", "PIN8", "PIN1", "PIN2", "PIN3", "PIN4", "AUX2", "THC1", "THC2", "PRESS"],
        "HSI0": 12,
        "HSI1": 7,
        "HSTime": 1,
        "HMStartEpoch": 1405087717.12159,
        "HMStartTime": "2014-07-11T16:08:37.121591+0200",
        "Info": {"DetectorInfo":"VxH:detbox=14952.235x0.000x1.000,dettab=-62.000x-245.000",
                 "ExperimentInfo":"0",
                 "MachineInfo": "Ie=183.27mA,u35u=100.000mm/0.001mm,u21m=100.000mm/0.000mm,u21d=100.000mm/-0.000mm",
                 "MirrorInfo": "rz=-3.600mrad,ty=0.455mm,ty1=2.075mm,ty2=-1.165mm,tz1=-0.030mm,tz2=-0.090mm,mir2rz=2.000mrad",
                 "OpticsInfo": "egy=12460.0eV,theta=9.132deg,hgt=11.7mm,tilt=4.440deg,tra=1.963mm",
                 "ProposalInfo": 0,
                 "StationInfo": "ID02"
                 }
        }

"""
    TO_SKIP = ("entry", "hdf5_filename", "plugin_name")

    def __init__(self):
        Plugin.__init__(self)
        self.cycle = None
        self.c216 = None
        self.nxs = None
        self.hdf5_filename = None
        self.entry = None
        self.instrument = None
        self.group = None
        self.tfg_grp = None
        self.mcs_grp = None
        self.input2 = {}
        if "TANGO_HOST" not in os.environ:
            raise RuntimeError("No TANGO_HOST defined")

    def setup(self, kwargs=None):
        """
        see class documentation
        """
        Plugin.setup(self, kwargs)
        if "HS32F10" in self.input:
            self.input2.update(preproc(**self.input))
        else:
            self.input2.update(self.input)
        # for debugging
        self.input["input2"] = self.input2

        self.c216 = self.input2.get("c216", "id02/c216/0")
        self.cycle = self.input2.get("cycle", 1)
        if "hdf5_filename" not in self.input2:
            self.log_error("hdf5_filename not in input")
        self.hdf5_filename = self.input2.get("hdf5_filename")
        self.entry = self.input2.get("entry", "entry")
        self.instrument = self.input2.get("instrument", "ID02")

    def process(self):
        self.create_hdf5()
        self.read_c216()

    def create_hdf5(self):
        """
        Create a HDF5 file and data-structure
        """
        try:
            self.nxs = Nexus(self.hdf5_filename, mode="a", creator="dahu")
        except IOError as error:
            self.log_warning("Unable to open %s: %s. Removing file and starting from scratch" % (self.hdf5_filename, error))
            os.unlink(self.hdf5_filename)
            self.nxs = Nexus(self.hdf5_filename, mode="w", creator="dahu")

        entry = self.nxs.new_entry(self.entry,
                                   program_name=self.input.get("plugin_name", fully_qualified_name(self.__class__)),
                                   title="C216 metadata collection")
        self.entry = entry.name
        entry["program_name"].attrs["version"] = __version__

        # configuration
        config_grp = self.nxs.new_class(entry, "configuration", "NXnote")
        config_grp["type"] = "text/json"
        config_grp["data"] = json.dumps(self.input, indent=2, separators=(",\r\n", ": "))

        # Instrument
        instrument_grp = self.nxs.new_instrument(entry=entry, instrument_name=self.instrument)
        instrument_grp["name"] = "TruSAXS"
        self.instrument = instrument_grp.name

        # TimeFrameGenerator
        self.tfg_grp = self.nxs.new_class(instrument_grp, "TFG", "NXcollection")
        self.tfg_grp["device"] = str(self.c216)

        # MultiCounterScaler
        self.mcs_grp = self.nxs.new_class(instrument_grp, "MCS", "NXcollection")
        self.mcs_grp["device"] = str(self.c216)

        # Static metadata
        self.info_grp = self.nxs.h5.require_group(posixpath.join(self.instrument, "parameters"))
        self.info_grp.attrs["NX_class"] = "NXcollection"

        for field, value in self.input2.get("Info", {}).items():
            if field not in self.TO_SKIP and not isinstance(value, dict):
                try:
                    value.encode("ascii")
                except UnicodeEncodeError:
                    self.log_warning("Unicode Error in field %s: %s, skipping" % (field, value))
                except AttributeError as err:
                    self.log_warning("Attribute Error %s \n in field %s: %s, forcing to string." % (err, field, value))
                    self.info_grp[field] = str(value)
                else:
                    self.info_grp[field] = str(value)

        # Factor
        HS32F = self.input2.get("HS32F")
        if HS32F is not None:
            self.mcs_grp.create_dataset("HS32F", data=HS32F).attrs.__setitem__("interpretation", "spectrum")
        # Zero
        HS32Z = self.input2.get("HS32Z")
        if HS32Z is not None:
            self.mcs_grp.create_dataset("HS32Z", data=HS32Z).attrs.__setitem__("interpretation", "spectrum")
        # Name
        HS32N = self.input2.get("HS32N")
        if HS32N is not None:
            self.mcs_grp.create_dataset("HS32N", data=[i.encode() for i in HS32N]).attrs.__setitem__("interpretation", "spectrum")
        # Mode
        HS32M = self.input2.get("HS32M")
        if HS32M is not None:
            self.mcs_grp.create_dataset("HS32M", data=HS32M).attrs.__setitem__("interpretation", "spectrum")

        if HS32N and HS32Z and HS32F:
            self.mcs_grp.require_group("interpreted")

    def read_c216(self):
        """
        Manage the metadata coming from C216 Time Frame Generator
        """
        c216ds = PyTango.DeviceProxy(str(self.c216))
        if c216ds.CompStatus("Tango::RUNNING") == "Tango::ON":
            msg = "C216 is running while reading counters ... possible issue"
            self._logging.append(msg)
            logger.warning(msg)

        raw_status = c216ds.GetCompleteConfigAndStatus()
        status = {"TFU_MODE": raw_status[0],
                  "TFU_FRAME_MODE": raw_status[1],
                  "START_MODE": raw_status[2],
                  "FRAMES": raw_status[14] // 2,
                  "CYCLES": raw_status[17],
                  "CHANNELS": raw_status[19],
                  "ACQ_BANK": raw_status[65],
                  "TFU_CONFIG": raw_status[:28],
                  "DIGI_IN": raw_status[28:44],
                  "ANA_IN": raw_status[44:60]}

        cycles = status["CYCLES"]
        frames = status["FRAMES"]
        tfg = c216ds.GetTfuFrameTimes()
        self.tfg_grp["frame_time"] = tfg
        self.tfg_grp["frame_time"].attrs["interpretation"] = "scalar"
        self.tfg_grp["cycles"] = cycles
        self.tfg_grp["frames"] = frames
        # handle the case of multiple cycles: n*more frames, exposure time always the same
        tmp = numpy.outer(numpy.ones(cycles), tfg).ravel()
        live_time = tmp[1::2]
        self.tfg_grp["live_time"] = live_time
        self.tfg_grp["live_time"].attrs["interpretation"] = "scalar"
        self.tfg_grp["dead_time"] = tmp[0::2]
        self.tfg_grp["dead_time"].attrs["interpretation"] = "scalar"
        self.tfg_grp["delta_time"] = tmp.cumsum()[0::2]
        self.tfg_grp["delta_time"].attrs["interpretation"] = "scalar"
        for key in ["HMStartTime", "HMStartEpoch"]:
            if key in self.input2:
                if type(self.input2[key]) in StringTypes:
                    self.tfg_grp[key] = str(self.input2[key])
                else:
                    self.tfg_grp[key] = self.input2[key]
                self.tfg_grp[key].attrs["interpretation"] = "scalar"

        # raw scalers:
        raw_scalers = c216ds.ReadScalersForNLiveFrames([0, frames - 1])
        raw_scalers.shape = frames, -1
        counters = raw_scalers.shape[1]
        self.mcs_grp["HS32C"] = raw_scalers
        if "HSTime" in self.input2:
            pin = int(self.input2["HSTime"])
            if pin > counters:
                self.log_error("invalid pin number %s" % pin)
            self.mcs_grp["HSTime"] = pin
            self.mcs_grp["HSTime"].attrs["interpretation"] = "scalar"
            self.mcs_grp["HSTime"].attrs["counter"] = "1-based pin number"
            pin -= 1  # 1 based pin number
            time_counter = raw_scalers[:, pin]
            if "HS32F" in self.mcs_grp:
                factor = self.mcs_grp["HS32F"][pin]
            else:
                self.log_warning("No factors provided for time measurement: defaulting to 1e-6")
                factor = 1e-6
            measured_time = time_counter * factor
            self.mcs_grp["ExposureTime"] = measured_time
            self.mcs_grp["ExposureTime"].attrs["interpretation"] = "scalar"
        else:
            self.log_error("No HSTime pin number, using TFG time")
            measured_time = tfg[1::2]

        if ("HS32F" in self.mcs_grp) and ("HS32Z" in self.mcs_grp):
            #             if "interpreted" in self.mcs_grp:
            modes = numpy.zeros(counters, dtype=numpy.int32)
            if "HS32M" in self.mcs_grp:
                raw_modes = numpy.array(self.mcs_grp["HS32M"])
                modes[:raw_modes.size] = raw_modes
            # else: mode=0
            values = numpy.zeros((frames, counters), dtype=numpy.float32)
            exptime = numpy.outer(tfg[1::2], numpy.ones(counters))
            zero = numpy.outer(numpy.ones(frames), numpy.array(self.mcs_grp["HS32Z"]))
            factor = numpy.outer(numpy.ones(frames), numpy.array(self.mcs_grp["HS32F"]))
            values_int = (raw_scalers - zero * exptime) * factor
            values_avg = (raw_scalers / exptime - zero) * factor
            mask = (numpy.outer(numpy.ones(frames), modes)).astype(bool)
            nmask = numpy.logical_not(mask)
            values[mask] = values_avg[mask]
            values[nmask] = values_int[nmask]
            self.mcs_grp["HS32V"] = values.astype(numpy.float32)
            self.mcs_grp["HS32V"].attrs["interpretation"] = "scalar"
            for i, name in enumerate(self.mcs_grp["HS32N"]):
                fullname = "interpreted/%s" % name
                self.mcs_grp[fullname] = values[:, i]
                self.mcs_grp[fullname].attrs["interpretation"] = "scalar"

            sot = self.input2.get("ShutterOpeningTime", 0.0)
            sct = self.input2.get("ShutterClosingTime", 0.0)
            for name, value in (("ShutterOpeningTime", sot),
                                ("ShutterClosingTime", sct)):
                        self.mcs_grp[name] = value
                        self.mcs_grp[name].attrs["interpretation"] = "scalar"
            correction_time = (measured_time - sot + sct) / (measured_time - sot)

            for I in ("HSI0", "HSI1"):
                if I in self.input2:
                    dest = "Intensity" + I[-1]
                    pin = int(self.input2[I])
                    if pin > counters:
                        self.log_error("invalid pin number %s" % pin)
                    self.mcs_grp[I] = pin
                    self.mcs_grp[I].attrs["interpretation"] = "scalar"
                    self.mcs_grp[I].attrs["counter"] = "1-based pin number"
                    pin -= 1  # 1 based pin number got 0 based.
                    counter = values[:, pin]
#                     factor = self.mcs_grp["HS32F"][pin]
#                     zero = self.mcs_grp["HS32Z"][pin]
#                     measured = (counter - measured_time * zero) * factor
                    I_factor = float(self.input2.get(I + "Factor", 1.0))
                    self.mcs_grp[I + "Factor"] = I_factor
                    self.mcs_grp[I + "Factor"].attrs["interpretation"] = "scalar"
                    measured = counter * I_factor
                    self.mcs_grp[dest + "UnCor"] = measured
                    self.mcs_grp[dest + "UnCor"].attrs["interpretation"] = "scalar"
                    self.mcs_grp[dest + "ShutCor"] = measured * correction_time
                    self.mcs_grp[dest + "ShutCor"].attrs["interpretation"] = "scalar"
        else:
            self.log_error("Not factor/zero to calculate I0/I1")

    def teardown(self):
        self.output["c216_filename"] = self.hdf5_filename
        if self.group:
            self.output["c216_path"] = self.group.name
            # self.group.parent["end_time"] = str(get_isotime())
        if self.nxs:
            self.nxs.close()
        Plugin.teardown(self)
