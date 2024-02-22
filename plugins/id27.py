"""Online data analysis for ID27
"""

import os
import sys
import shutil
import shlex
import collections
import io
import logging
import re
import glob
import subprocess
import json
from threading import Semaphore
import fabio
import pyFAI
from dahu.plugin import Plugin
from dahu.factory import register
lock = Semaphore()
logger = logging.getLogger("id27")

RAW = "RAW_DATA"
PROCESSED = "PROCESSED_DATA"
WORKDIR = "/scratch/shared"
PREFIX = os.path.dirname(sys.executable) #where there scripts are

try:
    from cryio import crysalis
except ImportError:
    crysalis = None

try:
    from pyicat_plus.client.main import IcatClient
except ImportError:
    print("iCat connection will not work")
    IcatClient = None


def crysalis_config(calibration_path, calibration_name, number_of_frames,
                    omega_start, omega_step, center, distance, wave_length, exposure_time):
    """
    :param calibration_path: typically "/users/opid27/file_conversion/"
    :param calibration_name: typically "scan0001"
    :param number_of_frames: integrer
    :param omega_start: in deg ?
    :param omega_step: in deg ?
    :param center: in pixel ?
    :param distance: in mm ?
    :param wave_length: in angstrom
    :param exposure_time: in seconds
    :return 2 dicts, one containing the path crysalis_files, the second with scans
    """
    calibration_name = calibration_name or ""
    calibration_path = calibration_path or ""
    par_file = os.path.join(calibration_path, calibration_name + '.par')
    if not os.path.exists(par_file):
        par_file = "/users/opid27/file_conversion/scan0001.par"
    crysalis_files = {
        'par_file': par_file,
        'set_file': '/users/opid27/file_conversion/scan0001.set',
        'ccd_file': '/users/opid27/file_conversion/scan0001.ccd'}

    scans = collections.OrderedDict()
    scans[0] = [{'count': number_of_frames,
                 'omega': 0,
                 'omega_start': omega_start,
                 'omega_end': omega_start + number_of_frames * omega_step,
                 'pixel_size': 0.075,
                 'omega_runs': None,
                 'theta': 0,
                 'kappa': 0,
                 'phi': 0,
                 'domega': omega_step,
                 'dtheta': 0,
                 'dkappa': 0,
                 'dphi': 0,
                 'center_x': center[0],
                 'center_y': center[1],
                 'alpha': 50,
                 'dist': distance,
                 'l1': wave_length,
                 'l2': wave_length,
                 'l12': wave_length,
                 'b': wave_length,
                 'mono': 0.99,
                 'monotype': 'SYNCHROTRON',
                 'chip': [1024, 1024],
                 'Exposure_time': exposure_time,
            }]

    return crysalis_files, scans


def copy_set_ccd(crysalis_files, crysalis_dir, basename):
    "copy .set and .ccd file"
    shutil.copy(crysalis_files['set_file'], os.path.join(crysalis_dir, basename + '.set'))
    shutil.copy(crysalis_files['ccd_file'], os.path.join(crysalis_dir, basename + '.ccd'))


def unpack_processed(completed_process):
    """Convert a CompletedProcess object in to something which is serialisable
    
    :param completed_process: Return of subprocess.run
    :return: dict with the same content
    """
    res = {}
    for key in dir(completed_process):

        if key.startswith("_"):
            continue
        val = completed_process.__getattribute__(key)
        if callable(val):
            continue
        if "decode" in dir(val):
            val = val.decode()
        res[key] = val
    return res


def createCrysalis(scans, crysalis_dir, basename):
    """Unused"""
    runHeader = crysalis.RunHeader(basename.encode(), crysalis_dir.encode(), 1)
    runname = os.path.join(crysalis_dir, basename)
    runFile = []

    for omega_run in scans[0]:
        dscr = crysalis.RunDscr(0)
        dscr.axis = crysalis.SCAN_AXIS['OMEGA']
        dscr.kappa = omega_run['kappa']
        dscr.omegaphi = 0
        dscr.start = omega_run['omega_start']
        dscr.end = omega_run['omega_end']
        dscr.width = omega_run['domega']
        dscr.todo = dscr.done = omega_run['count']
        dscr.exposure = 1
        runFile.append(dscr)

    crysalis.saveRun(runname, runHeader, runFile)
    crysalis.saveCrysalisExpSettings(crysalis_dir)


def create_par_file(crysalis_files, crysalis_dir, basename):
    """Create a new par-file by copying a reference one and 
       changing the ....
    """    
    ref_par = crysalis_files['par_file']
    new_par = os.path.join(crysalis_dir, basename + '.par')
    if os.path.exists(new_par):
        "make a backup"
        os.rename(new_par, new_par+".bak")
    with io.open(new_par, 'w', encoding='iso-8859-1') as new_file:
        with io.open(ref_par, 'r', encoding='iso-8859-1') as old_file:
            for line in old_file:
                if line.startswith("FILE CHIP"):
                    new_file.write('FILE CHIP "' + basename + '.ccd" \n')
                else:
                    new_file.write(line)


def create_rsync_file(filename, folder="esp"):
    """create 2 files: 
        * .source with the LIMA-HDF5 file 
        * a bash script to synchronise back the reduced crysalis dataset
    
    :param filename: name of the LIMA-HDF5 file with the images
    :param folder: name of the folder to create. If None, just return the path of the script
    :return: the path of the crysalis directory  
    """
    splitted = filename.split("/")
    if RAW in splitted:
        splitted[splitted.index(RAW)] = PROCESSED
        destname =  "/".join(splitted)
        ddir = os.path.dirname(destname)
        if not os.path.isdir(ddir):
            with lock:
                if not os.path.isdir(ddir):
                    os.makedirs(ddir)
    else:
        destname = filename
    dest_dir = "/".join(splitted[3:-1])  # strip /data/visitor ... filename.h5
    dest_dir = os.path.join(WORKDIR, dest_dir)
    script = os.path.join(dest_dir, "sync")
    if folder is None:
        return script

    crysalis_dir = os.path.join(dest_dir, folder)
    if not os.path.exists(crysalis_dir):
        with lock:
            if not os.path.exists(crysalis_dir):
                os.makedirs(crysalis_dir)
    with open(os.path.join(dest_dir, ".source"), "a") as source:
        source.write(filename + os.linesep)

    if os.path.exists(script):
        with open(script, "a") as source:
            source.write(os.linesep.join([f'rsync -avx {shlex.quote(os.path.join(dest_dir, folder))} {shlex.quote(os.path.dirname(destname))}', ""]))
    else:
        with open(script, "w") as source:
            source.write(os.linesep.join(['#!/bin/sh', f'rsync -avx {shlex.quote(os.path.join(dest_dir, folder))} {shlex.quote(os.path.dirname(destname))}', '']))
        os.chmod(script, 0o755)

    return crysalis_dir


def crysalis_conversion(wave_length=None, distance=None,
                        center=(None, None),
                        omega_start=None, omega_step=None,
                        exposure_time=None,
                        number_of_points=None,
                        file_source_path=None,
                        scan_name=None,
                        calibration_path="", calibration_name="",
                        update_mask=False, update_par=False,
                        **kwargs):
    """Run the `eiger2crysalis` script synchronously
     
    """

    assert crysalis, "cryio is not installed"

    logger.info('start data convertion with parameters:\n' +
                f'wave_length={wave_length}\n' +
                f'distance={distance}\n' +
                f'center={center}\n' +
                f'omega_start={omega_start}\n' +
                f'omega_step={omega_step}\n' +
                f'file_source_path={file_source_path}\n' +
                f'scan_name={scan_name}')

    script_name = os.path.join(PREFIX, 'eiger2crysalis')

    crysalis_dir = create_rsync_file(file_source_path)

    parameters = [script_name,
                  "-w", f'{wave_length}',
                  "-d", f'{distance}',
                  "-b", f'{center[0]}', f'{center[1]}',
                  f"--omega={omega_start}+{omega_step}*index",
                  file_source_path,
                  "-o", os.path.join(crysalis_dir, f'{scan_name}_1_''{index}.esperanto')]
    if update_mask:
        parameters.append("--calc-mask")
    logger.info('starts with parameters: %s', parameters)

    res = subprocess.run(parameters, capture_output=True, check=False)
    crysalis_files, scans = crysalis_config(calibration_path, calibration_name, number_of_points, omega_start, omega_step, center, distance, wave_length, exposure_time)

    if not update_mask:
        copy_set_ccd(crysalis_files, crysalis_dir, scan_name)
    #createCrysalis(scans, crysalis_dir, scan_name) # .run file already implemented in eiger2crysalis
    if not update_par:
        create_par_file(crysalis_files, crysalis_dir, scan_name)

    return res


def crysalis_conversion_fscannd(wave_length=None,
                                distance=None,
                                center=(None, None),
                                omega_start=None,
                                omega_step=None,
                                exposure_time=None,
                                npoints=None,
                                motor_mode=None,
                                dirname=None,
                                scan_name=None,
                                calibration_path="",
                                calibration_name="",
                                **kwargs):
    assert crysalis, "cryio is not installed"
    script_name = os.path.join(PREFIX, 'eiger2crysalis')
    pattern = re.compile('eiger_([0-9]+).h5')
    filenames_to_convert = glob.glob(f'{dirname}/eiger_????.h5')
    results = {}

    for filepath in sorted(filenames_to_convert):

        filename = os.path.basename(filepath)
        g = pattern.match(filename)
        if g:
            number = int(g.group(1))
            crysalis_folder_name = f'esp_{number:04d}'
            crysalis_dir = create_rsync_file(filepath, crysalis_folder_name)

            if motor_mode == "ZIGZAG" and (number % 2):
                revert_omega_start = omega_start + (omega_step * npoints)
                omega_par = f"--omega={revert_omega_start}-{omega_step}*index"
            else:
                omega_par = f"--omega={omega_start}+{omega_step}*index"

            parameters = [script_name,
                          "-w", f'{wave_length}',
                          "-d", f'{distance}',
                          "-b", f'{center[0]}', f'{center[1]}',
                          omega_par,
                          filepath,
                          "-o", os.path.join(crysalis_dir,
                                             f'{scan_name}_{number}_1_''{index}.esperanto')]
            # print('starts with parameters:',parameters)
            results[filename] = subprocess.run(parameters, capture_output=True, check=False)
            crysalis_files, scans = crysalis_config(calibration_path, calibration_name, npoints, omega_start, omega_step, center, distance, wave_length, exposure_time)

            crysalis_scan_name = scan_name + '_' + str(number)
            copy_set_ccd(crysalis_files, crysalis_dir, crysalis_scan_name)
            createCrysalis(scans, crysalis_dir, crysalis_scan_name)
            create_par_file(crysalis_files, crysalis_dir, crysalis_scan_name)
    return results


def fabio_conversion(file_path,
                     scan_number,
                     folder="xdi",
                     fabioimage="tifimage",
                     extension="tif",
                     export_icat=False):
    "Convert a set of eiger files to cbf or tiff"
    results = []
    filename = os.path.join(file_path, scan_number, 'eiger_????.h5')
    files = {f:fabio.open(f).nframes for f in sorted(glob.glob(filename))} #since python 3.7 dict are ordered !
    all_single = max(files.values())==1
    file_path = file_path.rstrip("/")

    splitted = file_path.split("/")
    dset_name = splitted[-1]
    if RAW in splitted:
        raw_pos = splitted.index(RAW)
        splitted[raw_pos] = PROCESSED
        splitted.append(scan_number)
        splitted.insert(0, "/")
        dest_dir = os.path.join(*splitted)    
        sample_name = splitted[raw_pos + 1]
    else:
        dest_dir = os.path.join(file_path, scan_number)
        sample_name = "unknown sample"

    if len(files) == 0:
        raise RuntimeError(f"No such file {filename}")

    cnt = 0
    for idx_file, filename in enumerate(files):
        img_data_fabio = fabio.open(filename)
        basename = folder if (len(files) == 1 or all_single) else f"{folder}_{idx_file+1:04d}"

        dest_dir2 = os.path.join(dest_dir, basename)
        if not os.path.exists(dest_dir2):
            with lock:
                if not os.path.exists(dest_dir2):
                    os.makedirs(dest_dir2)
        for i in range(img_data_fabio.nframes):
            cnt += 1
            frame = img_data_fabio.getframe(i)
            conv = frame.convert(fabioimage)
            data = conv.data.astype('int32')
            data[data < 0] = 0
            conv.data = data
            output = os.path.join(dest_dir2, f"{dset_name}_{cnt if all_single else i+1:04d}.{extension}")
            conv.write(output)
            results.append(output)
    if export_icat:
        metadata = {"definition": "conversion",
                    "Sample_name": sample_name}
        try:
            send_icat(file_path, dest_dir2, metadata=metadata)
        except Exception as err:
            import traceback
            print(f"Error {type(err)}: {err}")
            traceback.print_exc(err, file=sys.stdout)
    return results

##########################
# Start the server part
##########################


@register
class CrysalisConversion(Plugin):
    """
    This is the basic plugin of cysalis conversion
       
    Typical JSON file:

{"wave_length": 0.3738,
 "distance": 196.44,
 "center": [1560.9227, 1682.7944],
 "omega_start":-32,
 "omega_step": 0.5,
 "exposure_time":0.1,
 "number_of_points":128,
 "file_source_path": "/data/id27/inhouse/blc13357/id27/Vanadinite_Weck/Vanadinite_Weck_0001/scan0001/eiger_0000.h5",
 "scan_name": "scan0001",
 "calibration_path": "/data/id27/inhouse/blc13357/id27/Vanadinite_Weck/Vanadinite_Weck_0001/scan0001/esp",
 "calibration_name":"scan0001",
 "plugin_name": "id27.CrysalisConversion",
 "update_mask": 0,
 "update_par": 0
}    
    """

    def process(self):
        Plugin.process(self)
        if not self.input:
            logger.error("input is empty")

        result = crysalis_conversion(**self.input)
        self.output["results"] = unpack_processed(result)

        script = create_rsync_file(self.input["file_source_path"], None)
        sync_results = subprocess.run([script], capture_output=True, check=False)
        self.output["sync"] = unpack_processed(sync_results)
        #TODO: send to icat


@register
class CrysalisConversionFscannd(Plugin):
    """
    This is the plugin of cysalis conversion for fscannd type scans

    Typical JSON file:
    {"wave_length": 1.54,
     "distance": 100,
     "center": [400,500],
     "omega_start":-90,
     "omega_step": 1,
     "exposure_time":0.1,
     "npoints":180,
     "motor_mode":"ZIGZAG",
     "dirname": "/data/id27/inhouse/sample",
     "scan_name": "scan",
     "calibration_path":"/data/id27/inhouse/sample/vanadinite",
     "calibration_name":"vanadinite"
    }
    
    """

    def process(self):
        Plugin.process(self)
        if not self.input:
            logger.error("input is empty")
        result = crysalis_conversion_fscannd(**self.input)
        self.output["results"] = {k: unpack_processed(v) for k, v in result.items()}
        random_filename = glob.glob(f'{self.input["dirname"]}/eiger*.h5')[0]
        script = create_rsync_file(random_filename, None)
        sync_results = subprocess.run([script], capture_output=True, check=False)
        self.output["sync"] = unpack_processed(sync_results)
        #TODO: send to icat


@register
class XdiConversion(Plugin):
    """
    This is the plugin to convert an HDF5 to a stack of TIFF files
       
    Typical JSON file:
    {"file_path": "/data/id27/inhouse/some/file.h5",
     "scan_number": "0001"
    }
    
    """

    def process(self):
        Plugin.process(self)
        if not self.input:
            logger.error("input is empty")

        file_path = self.input["file_path"]
        scan_number = self.input["scan_number"]
        results = fabio_conversion(file_path,
                                   scan_number,
                                   folder="xdi",
                                   fabioimage="tifimage",
                                   extension="tif",
                                   export_icat=False)
        self.output["output"] = results


@register
class Average(Plugin):
    """
    This is the plugin to average out an HDF5 stack of images
       
    Typical JSON file:
    {"file_path": "/data/id27/inhouse/some/directory",
     "scan_number": "scan0001"
    }
    # optional: "output_directory": "/path/to/some/output"
    """

    def process(self):
        Plugin.process(self)

        results = {}
        outputs = []
        basename = 'sum.h5'
        prefix = basename.split('.')[0]

        if not self.input:
            logger.error("input is empty")
        output_directory = self.input.get("output_directory")
        file_path = self.input["file_path"]
        scan_number = self.input["scan_number"]
        filename = os.path.join(file_path, scan_number, 'eiger_????.h5')
        filenames = sorted(glob.glob(filename))

        if len(filenames) == 0:
            raise RuntimeError(f"File does not exist {filename}")
        sample_name = "undefined sample"
        for idx_h5, filename in enumerate(filenames):
            if output_directory:
                dest_dir = output_directory
                sample, dataset = file_path.strip().strip("/").split("/")[-2:]
                tmpname = basename if len(filenames) == 1 else f'{prefix}_{idx_h5+1:04d}.h5'
                basename = "_".join((sample, dataset, scan_number, tmpname))
            else:
                tmpname = prefix if len(filenames) == 1 else f'{prefix}_{idx_h5+1:04d}'
                splitted = file_path.split("/")
                if RAW in splitted:
                    raw_pos = splitted.index(RAW)
                    splitted[raw_pos] = PROCESSED
                    splitted.append(scan_number)
                    splitted.append(tmpname)
                    splitted.insert(0, "/")
                    dest_dir = os.path.join(*splitted)
                    sample_name = splitted[raw_pos+1]
                else:
                    dest_dir = os.path.join(file_path, scan_number, tmpname)
            if not os.path.exists(dest_dir):
                with lock:
                    if not os.path.exists(dest_dir):
                        os.makedirs(dest_dir)

            output = os.path.join(dest_dir, basename)
            command = [os.path.join(PREFIX, 'pyFAI-average'), 
                       '-m', 'sum', '-F', 'lima', '-o', output, filename]
            results[filename] = unpack_processed(subprocess.run(command, capture_output=True, check=False))
            outputs.append(output)
        self.output["output_filename"] = outputs
        self.output["conversion"] = results
        #send to icat
        metadata = {"definition": "sum",
                    "Sample_name": sample_name}
        try:
            send_icat(file_path, dest_dir, metadata=metadata)
        except Exception as err:
            import traceback
            print(f"Error {type(err)}: {err}")
            traceback.print_exc(err, file=sys.stdout)


@register
class XdsConversion(Plugin):
    """
    This is the plugin to convert an HDF5 to a stack of CBF files for XDS
       
    Typical JSON file:
    {"file_path": "/data/id27/inhouse/some/file.h5",
     "scan_number": "0001"
    }
    
    """

    def process(self):
        Plugin.process(self)
        if not self.input:
            logger.error("input is empty")

        file_path = self.input["file_path"]
        scan_number = self.input["scan_number"]
        results = fabio_conversion(file_path,
                                   scan_number,
                                   folder="cbf",
                                   fabioimage="cbfimage",
                                   extension="cbf",
                                   export_icat=False)
        self.output["output"] = results

@register
class DiffMap(Plugin):
    """
    This plugin performs the azimuthal integration of a set of data and populate a diffraction map

    Typical JSON input structure:
    { "plugin_name": "id27.diffmap",
      "ponifile": "/tmp/geometry.poni",
      "maskfile": "/tmp/mask.msk",
      "unit": "2th_deg",
      "npt": 2000,
      "file_path": "/data/id27/inhouse/some/path",
      "scan_number": "scan_0001",
      "slow_scan": 10,
      "fast_scan": 10
    }
    """
    def process(self):
        Plugin.process(self)
        if not self.input:
            logger.error("input is empty")

        file_path = self.input["file_path"]
        scan_number = self.input["scan_number"]
        dest_dir = os.path.join(file_path.replace(RAW, PROCESSED), scan_number)
        filename = os.path.join(file_path, scan_number, 'eiger_????.h5')
        files = sorted(glob.glob(filename))
        if not os.path.isdir(dest_dir):
            with lock:
                if not os.path.isdir(dest_dir):
                    os.makedirs(dest_dir)
        config = os.path.join(dest_dir, "diff_map.json")
        dest = os.path.join(dest_dir, "diff_map.h5")
        results = {}
        param = {}
        ai = pyFAI.load(self.input.get("ponifile")).get_config()
        if "maskfile" in self.input:
            ai["do_mask"] = True
            ai["mask_file"] = self.input["maskfile"]
        # some constants hardcoded for the beamline:
        ai["do_polarization"] = True
        ai["polarization_factor"] = 0.99
        ai["do_solid_angle"] = True
        ai["error_model"] = "poisson"
        ai["application"] = "pyfai-integrate"
        ai["version"] = 3
        ai["method"] = ["bbox", "csr", "opencl"]
        ai["opencl_device"] = "gpu"
        ai["nbpt_rad"] = self.input.get("npt", 1)
        ai["nbpt_azim"] = 1
        ai["do_2D"] = False
        ai["unit"] = self.input.get("unit", "q_nm^-1")
        param["ai"] = ai
        param["experiment_title"] = os.path.join(os.path.basename(file_path), scan_number)
        param["fast_motor_name"] = "fast"
        param["slow_motor_name"] = "slow"
        param["fast_motor_points"] = self.input.get("fast_scan", 1)
        param["slow_motor_points"] = self.input.get("slow_scan", 1)
        param["offset"] = 0
        param["output_file"] = dest
        param["input_data"] = [(i, None, None) for i in files]
        with open(config, "w") as w:
            w.write(json.dumps(param, indent=2))
        results["config"] = config
        command = [os.path.join(PREFIX, 'pyFAI-diffmap'), '--no-gui', '--config', config]
        results["processing"] = unpack_processed(subprocess.run(command, capture_output=True, check=False))
        self.output["output_filename"] = dest
        self.output["diffmap"] = results
        #send to icat
        if RAW in file_path:
            l = file_path.split("/")
            i = l.index(RAW)
            sample_name = l[i+1]
        metadata = {"definition": "diffmap",
                    "Sample_name": sample_name}
        try:
            send_icat(file_path, dest_dir, metadata=metadata)
        except Exception as err:
            import traceback
            print(f"Error {type(err)}: {err}")
            traceback.print_exc(err, file=sys.stdout)


def send_icat(raw_dir, processed_dir, beamline="id27", proposal="", dataset="", metadata=None):
    "Function that sends to icat the processed data"
    icat_client = IcatClient(metadata_urls=["bcu-mq-01.esrf.fr:61613", "bcu-mq-02.esrf.fr:61613"])
    metadata = metadata or {"definition": "dummy processing", "Sample_name": "unknown sample"}
    l = raw_dir.split("/")
    if not proposal:
        try:
            visitor_idx = l.index("visitor")
        except:
            proposal = "undefined"
        else:
            proposal = l[visitor_idx+1]
    if not dataset:
        dataset = l[-1]
    kwargs = {"beamline":beamline, 
              "proposal":proposal, 
              "dataset":dataset, 
              "path":processed_dir, 
              "metadata":metadata, 
              "raw":[raw_dir]}
    icat_client.store_processed_data(**kwargs)
    return kwargs
    
    
    
    
