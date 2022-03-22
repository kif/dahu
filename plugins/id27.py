import os
import shutil
import collections
import io
import logging
import re
import glob
from dahu.plugin import Plugin
from dahu.factory import register
import fabio
logger = logging.getLogger("id27")

try:
    from cryio import crysalis
except ImportError:
    crysalis = None

import subprocess
# import threading
# import multiprocessing

WORKDIR = "/scratch/shared"


def crysalis_config(calibration_path,calibration_name, number_of_frames,  omega_start, omega_step, center, distance, wave_length, exposure_time):
    crysalis_files = {
        'par_file': os.path.join(calibration_path,calibration_name+'.par'),
        'set_file': '/users/opid27/file_conversion/scan0001.set',
        'ccd_file': '/users/opid27/file_conversion/scan0001.ccd'}


    scans = collections.OrderedDict()
    scans[0] = [{
            'count': number_of_frames,
            'omega': 0,
            'omega_start': omega_start,
            'omega_end': omega_start+number_of_frames*omega_step,
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
            'chip': [1024,1024],
            'Exposure_time': exposure_time,
            }]

    return crysalis_files, scans

def copy_set_ccd(crysalis_files, crysalis_dir, basename):

    shutil.copy(crysalis_files['set_file'], os.path.join(crysalis_dir,basename+'.set'))
    shutil.copy(crysalis_files['ccd_file'], os.path.join(crysalis_dir,basename+'.ccd'))

def unpack_CompletedProcess(cp):
    """Convert a CompletedProcess object in to something which is serialisable
    
    :param cp: Return of subprocess.run
    :return: dict with the same content
    """
    return {k: cp.__getattribute__(k) for k in dir(cp) if not (k.startswith("_") or callable(cp.__getattribute__(k)))}
    

def createCrysalis(scans, crysalis_dir, basename):
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

    new_par = os.path.join(crysalis_dir,basename+'.par')

    with io.open(new_par, 'w',encoding='iso-8859-1') as new_file:
        with io.open(crysalis_files['par_file'], 'r',encoding='iso-8859-1') as old_file:
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
    dest_dir = "/".join(filename.split("/")[3:-1]) #strip /data/visitor ... filename.h5
    dest_dir = os.path.join(WORKDIR, dest_dir)
    script = os.path.join(dest_dir, "sync")
    if folder is None:
        return script

    crysalis_dir = os.path.join(dest_dir, folder)
    if not os.path.exists(crysalis_dir):
            os.makedirs(crysalis_dir)
    with open(os.path.join(dest_dir, ".source"), "a") as source:
        source.write(filename+os.linesep)
    
    if os.path.exists(script):
        with open(script, "a") as source:
            source.write(os.linesep.join([f'rsync -avx {os.path.join(dest_dir, folder)} {os.path.dirname(filename)}', ""]))
    else:
        with open(script, "w") as source:
            source.write(os.linesep.join(['#!/bin/sh', f'rsync -avx {os.path.join(dest_dir, folder)} {os.path.dirname(filename)}', '']))
        os.chmod(script, 0o755)

    return crysalis_dir


def crysalis_conversion(wave_length=None, distance=None,
                        center=(None, None),
                        omega_start=None,omega_step=None,
                        exposure_time=None,
                        number_of_points=None,
                        file_source_path=None,
                        scan_name=None,
                        calibration_path="", calibration_name="",
                        **kwargs):
    """Run the `eiger2crysalis` script synchronously
     
    """
    
    assert crysalis, "cryio is not installed"
    
    logger.info('start data convertion with parameters:\n'+
                f'wave_length={wave_length}\n'+
                f'distance={distance}\n'+
                f'center={center}\n'+
                f'omega_start={omega_start}\n'+
                f'omega_step={omega_step}\n'+
                f'file_source_path={file_source_path}\n'+
                f'scan_name={scan_name}')

    script_name = 'eiger2crysalis'

    crysalis_dir = create_rsync_file(file_source_path)
    
    parameters = [script_name,
                  "-w", f'{wave_length}',
                  "-d", f'{distance}',
                  "-b", f'{center[0]}', f'{center[1]}',
                  f"--omega={omega_start}+{omega_step}*index",
                  file_source_path,
                  "-o", os.path.join(crysalis_dir, f'{scan_name}_1_''{index}.esperanto')]
    logger.info('starts with parameters: %s',parameters)

    crysalis_files, scans = crysalis_config(calibration_path, calibration_name, number_of_points,  omega_start, omega_step, center, distance, wave_length, exposure_time)

    copy_set_ccd(crysalis_files,crysalis_dir, scan_name )
    createCrysalis(scans, crysalis_dir, scan_name)
    create_par_file(crysalis_files,crysalis_dir, scan_name)

    return subprocess.run(parameters, capture_output=True)


def crysalis_conversion_fscannd(wave_length=None,
                                distance=None,
                                center=(None,None),
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
    script_name = 'eiger2crysalis'
    pattern = re.compile('eiger_([0-9]+).h5')
    filenames_to_convert = glob.glob(f'{dirname}/eiger*.h5')
    results = {}
    if filenames_to_convert:
        crysalis
    for filepath in sorted(filenames_to_convert):
        
        filename = os.path.basename(filepath)
        g = pattern.match(filename)
        if g:
            number = int(g.group(1))
            crysalis_folder_name = f'esp_{number}'
            crysalis_dir = create_rsync_file(filepath, crysalis_folder_name)
            
            if motor_mode == "ZIGZAG" and (number % 2):
                revert_omega_start=omega_start+(omega_step*npoints)
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
            #print('starts with parameters:',parameters)
            results[filename] = str(subprocess.run(parameters, capture_output=True))
            crysalis_files, scans = crysalis_config(calibration_path, calibration_name, npoints,  omega_start, omega_step, center, distance, wave_length, exposure_time)
            
            crysalis_scan_name = scan_name + '_' + str(number)
            copy_set_ccd(crysalis_files,crysalis_dir, crysalis_scan_name )
            createCrysalis(scans, crysalis_dir,  crysalis_scan_name)
            create_par_file(crysalis_files,crysalis_dir, crysalis_scan_name)
    return results


def fabio_conversion(file_path,
                     scan_number,
                     folder="xdi",
                     fabioimage="tifimage",
                     extension="tif"):

        filename = os.path.join(file_path,scan_number,'eiger_0000.h5')
        img_data_fabio = fabio.open(filename)
        dest_dir = os.path.join(file_path,scan_number, folder)
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)
        results = []
        for i, frame in enumerate(img_data_fabio):
            conv = frame.convert(fabioimage)
            data = conv.data.astype('int32')
            data[data <0 ] = 0
            conv.data = data
            #print('xdi converison', i)
            output = os.path.join(dest_dir, f"frame_{i+1:04d}.{extension}")
            conv.write(output)
            results.append(output)
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
 "plugin_name": "id27.CrysalisConversion"
}    
    """
    
    def process(self):
        Plugin.process(self)
        if not self.input:
            logger.error("input is empty")

        result = crysalis_conversion(**self.input)
        self.output["results"] = unpack_CompletedProcess(result)   

        script = create_rsync_file(self.input["file_source_path"], None)
        sync_results = subprocess.run([script], capture_output=True)
        self.output["sync"] = unpack_CompletedProcess(sync_results)


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
        self.output["results"] = {k: unpack_CompletedProcess(v) for k,v in result.items()}
        random_filename = glob.glob(f'{self.input["dirname"]}/eiger*.h5')[0]        
        script = create_rsync_file(random_filename, None)
        sync_results = subprocess.run([script], capture_output=True)
        self.output["sync"] = unpack_CompletedProcess(sync_results)


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
                         extension="tif")
        self.output["output"] = results        

@register
class Average(Plugin):
    """
    This is the plugin to average out an HDF5 stack of images
       
    Typical JSON file:
    {"file_path": "/data/id27/inhouse/some/directory",
     "scan_number": "scan0001"
    }
    
    """
    
    def process(self):
        Plugin.process(self)
        if not self.input:
            logger.error("input is empty")
            
            

        file_path = self.input["file_path"]
        scan_number = self.input["scan_number"]
        filename = os.path.join(file_path, scan_number,'eiger_0000.h5')

        dest_dir = os.path.join(file_path, scan_number, 'sum')

        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)

        output = os.path.join(dest_dir,'sum.edf')
        command = ['pyFAI-average', '-m', 'sum', '-o', output, filename ]
        result = subprocess.run(command, capture_output=True)
        self.output["output_filename"] = output
        self.output["results"] = str(result)        

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
                         extension="cbf")
        self.output["output"] = results        
