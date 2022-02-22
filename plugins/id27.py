import os
import shutil
import collections
import io
import logging
from dahu.plugin import Plugin
from dahu.factory import register
import fabio
logger = logging.getLogger("id27")

from cryio import crysalis
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


def create_rsync_file(filename):
    """create 2 files: 
        * .source with the LIMA-HDF5 file 
        * a bash script to synchronise back the reduced crysalis dataset
    
    :param filename: name of the LIMA-HDF5 file with the images
    :return: the path of the crysalis directory  
    """
    dest_dir = "/".join(filename.split("/")[3:-1]) #strip /data/visitor ... filename.h5
    dest_dir = os.path.join(WORKDIR, dest_dir)
    crysalis_dir = os.path.join(dest_dir, 'esp')
    if not os.path.exists(crysalis_dir):
            os.makedirs(crysalis_dir)
    with open(os.path.join(dest_dir, ".source"), "w") as source:
        source.write(filename)
    with open(os.path.join(dest_dir, "sync"), "w") as source:
        source.write(f'#!/bin/sh\nrsync -avx {os.path.join(dest_dir, "esp")} {os.path.dirname(filename)}\n')
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

        return subprocess.run(parameters)

##########################
# Start the server part
##########################

@register
class CrysalisConversion(Plugin):
    """
    This is the basic plugin of cysalis conversion
       
    Typical JSON file:
    {"wave_length": 1.54,
     "distance": 100,
     "center": [400,500],
     "omega_start":-90,
     "omega_step": 1,
     "exposure_time":0.1,
     "number_of_points":180,
     "file_source_path": "/data/id27/inhouse/sample/example.h5",
     "scan_name": "scan",
     "calibration_path":"/data/id27/inhouse/sample/vanadinite",
     "calibration_name":"vanadinite"
    }
    
    """
    
    def process(self):
        Plugin.process(self)
        if self.input is None:
            logger.warning("input is None")
        result = crysalis_conversion(**self.input)
        self.output["results"] = result        
