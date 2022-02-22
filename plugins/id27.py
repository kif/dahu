import logging
from dahu.plugin import Plugin
from dahu.factory import register

logger = logging.getLogger("id15")

import os, re, glob, sys, shutil, collections, io
from cryio import crysalis
import subprocess
import threading
import multiprocessing

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


def create_par_file(crysalis_files,crysalis_dir, basename):

    new_par = os.path.join(crysalis_dir,basename+'.par')

    with io.open(new_par, 'w',encoding='iso-8859-1') as new_file:
        with io.open(crysalis_files['par_file'], 'r',encoding='iso-8859-1') as old_file:
            for line in old_file:
                if line.startswith("FILE CHIP"):
                    new_file.write('FILE CHIP "' + basename + '.ccd" \n')
                else:
                    new_file.write(line)



########################################################## Start the server part

def crysalis_conversion(wave_length=None,distance=None,
                            center=(None,None),
                            omega_start=None,omega_step=None,
                            exposure_time=None,
                            number_of_points=None,
                            file_source_path=None,
                            scan_name=None,
                            calibration_path="",calibration_name="", **kwargs):
        print('start data convertion with parameters:\n',
              f'wave_length={wave_length}\n'
              f'distance={distance}\n',
              f'center={center}\n',
              f'omega_start={omega_start}\n',
              f'omega_step={omega_step}\n',
              f'file_source_path={file_source_path}\n',
              f'scan_name={scan_name}')

        script_name = 'eiger2crysalis'

        dirname = os.path.split(file_source_path)
        parameters = [script_name,
                      "-w", f'{wave_length}',
                      "-d", f'{distance}',
                      "-b", f'{center[0]}', f'{center[1]}',
                      f"--omega={omega_start}+{omega_step}*index",
                      file_source_path,
                      "-o", os.path.join(dirname[0],'esp',f'{scan_name}_1_''{index}.esperanto')]
        print('starts with parameters:',parameters)

        crysalis_files, scans = crysalis_config(calibration_path, calibration_name, number_of_points,  omega_start, omega_step, center, distance, wave_length, exposure_time)
        crysalis_dir = os.path.join(dirname[0], 'esp')
        isExist = os.path.exists( crysalis_dir)
        if not isExist:
            os.makedirs(crysalis_dir)
        copy_set_ccd(crysalis_files,crysalis_dir, scan_name )
        createCrysalis(scans, crysalis_dir, scan_name)
        create_par_file(crysalis_files,crysalis_dir, scan_name)

        return subprocess.run(parameters)


@register
class CrysalisConversion(Plugin):
    def process(self):
        Plugin.process(self)
        if self.input is None:
            logger.warning("input is None")
        result = crysalis_conversion(**self.input)
        self.log_warning(str(result))
        
                