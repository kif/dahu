#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Data Analysis plugin for BM29: BioSaxs

* SubtractBuffer: Search for the equivalence of buffers, average them and subtract from sample signal.  
* SaxsAnalysis: Performs Guinier + Kratky + IFT, generates plots  
"""

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "27/05/2020"
__status__ = "development"
__version__ = "0.1.0"

import os
import json
from math import log, pi
from collections import namedtuple
from dahu.plugin import Plugin
from dahu.factory import register
from dahu.cache import DataCache
from dahu.utils import fully_qualified_name
import logging
logger = logging.getLogger("bm29.subtract")
import numpy
try:
    import numexpr
except ImportError:
    logger.error("Numexpr is not installed, falling back on numpy's implementations")
    numexpr = None
import h5py
import fabio
import pyFAI, pyFAI.azimuthalIntegrator
from pyFAI.method_registry import IntegrationMethod
import freesas, freesas.cormap
from freesas.autorg import auto_gpa, autoRg
from freesas.bift import BIFT
from scipy.optimize import minimize
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from .common import Sample, Ispyb, get_equivalent_frames, cmp_float, get_integrator, KeyCache, polarization_factor, method, Nexus, get_isotime

NexusJuice = namedtuple("NexusJuice", "filename h5path npt unit q I poni mask energy polarization method signal2d error2d buffer concentration")


class SubtractBuffer(Plugin):
    """Search for the equivalence of buffers, average them and subtract from sample signal.
    
        Typical JSON file:
    {
      "buffers_files" = ["buffer_001.h5", "buffer_002.h5"],
      "sample_file" = "sample.h5",
      "output_file" = "subtracted.h5"
      "fidelity" = 0.001,
      # TODO ... "wait_for" = ["job1", "job2"]
      
    } 
    """
    def __init__(self):
        Plugin.__init__(self)
        self.buffer_files = []
        self.sample_file = None
        self.nxs = None
        self.output_file = None
        self.npt = None
        self.poni = None
        self.mask = None
        self.energy = None
        self.sample_juice = None
        self.buffer_juices = []
        self.Rg = self.I0 = self.Dmax = self.Vc = self.mass = None
           
    def setup(self, kwargs=None):
        logger.debug("SubtractBuffer.setup")
        Plugin.setup(self, kwargs)
        self.sample_file = self.input.get("sample_file")
        if self.sample_file is None:
            self.log_error("No sample file provided", do_raise=True)
        self.output_file = self.input.get("output_file")
        if self.output_file is None:
            lst = list(os.path.splitext(self.sample_file))
            lst.insert(1, "-sub")
            self.output_file = "".join(lst)
            self.log_warning(f"No output file provided, using: {self.output_file}")
        self.buffer_files = [os.path.abspath(fn) for fn in self.input.get("buffer_files", [])
                             if os.path.exists(fn)]
            
    def teardown(self):
        Plugin.teardown(self)
        logger.debug("SubtractBuffer.teardown")
        # export the output file location
        self.output["output_file"] = self.output_file
        self.output["Rg"] = self.Rg
        self.output["I0"] = self.I0
        self.output["Dmax"] = self.Dmax
        self.output["Vc"] = self.Vc
        self.output["mass"] = self.mass
        if self.nxs is not None:
            self.nxs.close()
        self.sample_juice = None
        self.buffer_juices = []

    def process(self):
        Plugin.process(self)
        logger.debug("SubtractBuffer.process")
        self.sample_juice = self.read_nexus(self.sample_file)
        self.create_nexus()
        
    def validate_buffer(self, buffer_file):
        "Validate if a buffer is consitent with the sample, return some buffer_juice or None when unconsistent"
        buffer_juice = self.read_nexus(buffer_file)
        if self.sample_juice.npt != buffer_juice.npt:
            self.log_warning(f"Sample {buffer_file} differs in number of points, discarding")
            return
        if abs(self.sample_juice.q - buffer_juice.q).max() > 1e-6:
            self.log_warning(f"Sample {buffer_file} differs in q-position, discarding")
            return
        if self.sample_juice.poni != buffer_juice.poni:
            self.log_warning(f"Sample {buffer_file} differs in poni-file, discarding")
            return
        if self.sample_juice.mask != buffer_juice.mask:
            self.log_warning(f"Sample {buffer_file} differs in mask-file, discarding")
            return
        if self.sample_juice.polarization != buffer_juice.polarization:
            self.log_warning(f"Sample {buffer_file} differs in polarization factor, discarding")
            return
        if self.sample_juice.buffer != buffer_juice.buffer:
            self.log_warning(f"Sample {buffer_file} differs in buffer descsription, discarding")
            return
        if buffer_juice.concentration:
            self.log_warning(f"Sample {buffer_file} concentration not null, discarding")
            return
        return buffer_juice
    
    def create_nexus(self):
        nxs = self.nxs = Nexus(self.output_file, mode="w")
        entry_grp = nxs.new_entry("entry", self.input.get("plugin_name", "dahu"), 
                              title='BioSaxs buffer subtraction', 
                              force_time=get_isotime())
        nxs.h5.attrs["default"] = entry_grp.name
        
    # Configuration
        cfg_grp = nxs.new_class(entry_grp, "configuration", "NXnote")
        cfg_grp.create_dataset("data", data=json.dumps(self.input, indent=2, separators=(",\r\n", ":\t")))
        cfg_grp.create_dataset("format", data = "text/json")

    #Process 0: Measurement group
        input_grp = nxs.new_class(entry_grp, "0_measurement", "NXcollection")
        input_grp["sequence_index"] = 0
        rel_path = os.path.relpath(os.path.abspath(self.sample_file), os.path.dirname(os.path.abspath(self.output_file)))
        input_grp["sample"] = h5py.ExternalLink(rel_path, self.sample_juice.h5path)

        for idx, buffer_file in enumerate(self.buffer_files):
            buffer_juice = self.validate_buffer(buffer_file)
            if buffer_file is not None:
                rel_path = os.path.relpath(os.path.abspath(buffer_file), os.path.dirname(os.path.abspath(self.output_file)))
                input_grp["buffer_%i"%idx] = h5py.ExternalLink(rel_path, buffer_juice.h5path)     
                self.buffer_juices.append(buffer_juice)
        
    #Process 1: CorMap
        cormap_grp = nxs.new_class(entry_grp, "1_correlation_mapping", "NXprocess")
        cormap_grp["sequence_index"] = 1
        cormap_grp["program"] = "freesas.cormap"
        cormap_grp["version"] = freesas.version
        cormap_grp["date"] = get_isotime()
        cormap_data = nxs.new_class(cormap_grp, "results", "NXdata")
        cfg_grp = nxs.new_class(cormap_grp, "configuration", "NXcollection")
        
    #Cormap processing   
        nb_frames = len(self.buffer_juices)
        count = numpy.empty((nb_frames, nb_frames), dtype=numpy.uint16)
        proba = numpy.empty((nb_frames, nb_frames), dtype=numpy.float32)
        for i in range(nb_frames):
            proba[i, i] = 1.0
            count[i, i] = 0
            for j in range(i):
                res = freesas.cormap.gof(self.buffer_juices[i].I, self.buffer_juices[j].I)
                proba[i,j] = proba[j,i] = res.P
                count[i,j] = count[j,i] = res.c
        fidelity = self.input.get("fidelity", 0)
        cfg_grp["fidelity_abs"] = fidelity
        cfg_grp["fidelity_rel"] = fidelity
        tomerge = get_equivalent_frames(proba, fidelity, fidelity)

        cormap_data.attrs["signal"] = "probability"
        cormap_ds = cormap_data.create_dataset("probability", data=proba)
        cormap_ds.attrs["interpretation"] = "image"
        cormap_ds.attrs["long_name"] = "Probability to be the same"
        
        count_ds = cormap_data.create_dataset("count", data=count)
        count_ds.attrs["interpretation"] = "image"
        count_ds.attrs["long_name"] = "Longest sequence where curves do not cross each other"

        to_merge_ds = cormap_data.create_dataset("to_merge",data= numpy.arange(*tomerge, dtype=numpy.uint16))
        to_merge_ds.attrs["long_name"] = "Index of equivalent frames" 
        cormap_grp.attrs["default"] = cormap_data.name
        
    #Process 2: Image processing: subtraction with standard deviation
        average_grp = nxs.new_class(entry_grp, "2_buffer_subtraction", "NXprocess")
        average_grp["sequence_index"] = 2
        average_grp["program"] = fully_qualified_name(self.__class__)
        average_grp["version"] = __version__
        average_data = nxs.new_class(average_grp, "results", "NXdata")
        average_data.attrs["signal"] = "intensity_normed"
    # Stage 2 processing

        #Nota: formula probably wrong ! Take into account the number of readings !
        #TODO implement those math using numexpr:
        #  avg = Σdata / Σnorm
        #  var = Σdata / (Σnorm)²
        #  Σnorm = avg / var
        #  Σdata = avg² / var
        # The propagate the error based on the number of frames in each buffer with quadratic summation
        # 
        buffer_average = numpy.mean([self.buffer_juices[i].signal2d for i in range(*tomerge)], axis=0)
        buffer_variance = numpy.sum([(self.buffer_juices[i].error2d)**2 for i in range(*tomerge)], axis=0) / (tomerge[1] - tomerge[0])**2
        sub_average = self.sample_juice.signal2d - buffer_average
        sub_variance = self.sample_juice.error2d**2 + buffer_variance
        sub_std = numpy.sqrt(sub_variance)
        
        
        int_avg_ds =  average_data.create_dataset("intensity_normed", 
                                                  data=numpy.ascontiguousarray(sub_average, dtype=numpy.float32),
                                                  **cmp_float)
        int_avg_ds.attrs["interpretation"] = "image"
        int_avg_ds.attrs["formula"] = "sample_signal - mean(buffer_signal_i)"
        int_std_ds =  average_data.create_dataset("intensity_std", 
                                                   data=numpy.ascontiguousarray(sub_std, dtype=numpy.float32),
                                                   **cmp_float)
        int_std_ds.attrs["interpretation"] = "image"    
        int_std_ds.attrs["formula"] = "sqrt( sample_variance + (sum(buffer_variance)/n_buffer**2 ))"
        int_std_ds.attrs["method"] = "quadratic sum of sample error and buffer errors"
        average_grp.attrs["default"] = average_data.name
    # Process 3: Azimuthal integration of the subtracted image
        ai2_grp = nxs.new_class(entry_grp, "3_azimuthal_integration", "NXprocess")
        ai2_grp["sequence_index"] = 3
        ai2_grp["program"] = "pyFAI"
        ai2_grp["version"] = pyFAI.version
        ai2_grp["date"] = get_isotime()
        key_cache = KeyCache(self.sample_juice.npt, self.sample_juice.unit, self.sample_juice.poni, self.sample_juice.mask, self.sample_juice.energy)
        ai = get_integrator(key_cache)
        radial_unit, unit_name = str(key_cache.unit).split("_", 1)
        ai2_data = nxs.new_class(ai2_grp, "results", "NXdata")
        ai2_data.attrs["signal"] = "I"
        ai2_data.attrs["axes"] = radial_unit
        ai2_grp.attrs["default"] = ai2_data.name
        cfg_grp = nxs.new_class(ai2_grp, "configuration", "NXnote")
        cfg_grp.create_dataset("data", data=json.dumps(ai.get_config(), indent=2, separators=(",\r\n", ": ")))
        cfg_grp.create_dataset("format", data = "text/json")
        if os.path.exists(key_cache.poni):
            cfg_grp.create_dataset("file_name", data = key_cache.poni)
        pol_ds = cfg_grp.create_dataset("polarization_factor", data=polarization_factor)
        pol_ds.attrs["comment"] = "Between -1 and +1, 0 for circular"
        cfg_grp.create_dataset("integration_method", data=json.dumps(method.method._asdict()))

    # Stage 3 processing: azimuthal integration
        res2 = ai._integrate1d_ng(sub_average, key_cache.npt, 
                                  variance=sub_variance,
                                  polarization_factor=self.sample_juice.polarization,
                                  unit=key_cache.unit,
                                  safe=False,
                                  method=self.sample_juice.method)

        ai2_q_ds = ai2_data.create_dataset(radial_unit,
                                           data=numpy.ascontiguousarray(res2.radial, dtype=numpy.float32))
        ai2_q_ds.attrs["units"] = unit_name
        radius_unit = "nm" if "nm" in unit_name else "Å"
        ai2_q_ds.attrs["long_name"] = "Scattering vector q (nm⁻¹)"
        
        ai2_int_ds = ai2_data.create_dataset("I", data=numpy.ascontiguousarray(res2.intensity, dtype=numpy.float32))
        ai2_std_ds = ai2_data.create_dataset("errors", 
                                             data=numpy.ascontiguousarray(res2.sigma, dtype=numpy.float32))
        
        
        ai2_int_ds.attrs["interpretation"] ="spectrum"     
        ai2_int_ds.attrs["units"] = "arbitrary"
        ai2_int_ds.attrs["long_name"] = "Intensity (absolute, normalized on water)"
        #ai2_int_ds.attrs["uncertainties"] = "errors" #this does not work
        ai2_std_ds.attrs["interpretation"] ="spectrum"
        ai2_int_ds.attrs["units"] = "arbitrary"
    
    # Process 4: Guinier analysis
        guinier_grp = nxs.new_class(entry_grp, "4_Guinier_analysis", "NXprocess")
        guinier_grp["sequence_index"] = 4
        guinier_grp["program"] = "freesas.autorg"
        guinier_grp["version"] = freesas.version
        guinier_grp["date"] = get_isotime()
        guinier_autorg = nxs.new_class(guinier_grp, "autorg", "NXcollection")
        guinier_gpa = nxs.new_class(guinier_grp, "gpa", "NXcollection")
        guinier_data = nxs.new_class(guinier_grp, "results", "NXdata")
    # Stage4 processing: autorg and auto_gpa
        sasm = numpy.vstack((res2.radial, res2.intensity, res2.sigma)).T

        try:
            gpa = auto_gpa(sasm)
        except Exception as error:
            guinier_gpa["Failed"] = "%s: %s"%(error.__class__.__name__, error)
            gpa = None
        else:
            #"Rg sigma_Rg I0 sigma_I0 start_point end_point quality aggregated"
            guinier_gpa["Rg"] = gpa.Rg
            guinier_autorg["Rg"].attrs["unit"] = radius_unit
            guinier_gpa["Rg_error"] = gpa.sigma_Rg
            guinier_autorg["Rg_error"].attrs["unit"] = radius_unit
            guinier_gpa["I0"] = gpa.I0
            guinier_gpa["I0_error"] = gpa.sigma_I0
            guinier_gpa["start_point"] = gpa.start_point
            guinier_gpa["end_point"] = gpa.end_point
#             guinier_gpa["quality"] = autorg.quality
#             guinier_gpa["aggregated"] = autorg.aggregated
        try:
            autorg = autoRg(sasm)
        except Exception as err:
            guinier_autorg["Failed"] = "%s: %s"%(err.__class__.__name__, err)
            autorg = None
        else:
            #"Rg sigma_Rg I0 sigma_I0 start_point end_point quality aggregated"
            guinier_autorg["Rg"] = autorg.Rg
            guinier_autorg["Rg"].attrs["unit"] = radius_unit
            guinier_autorg["Rg_error"] = autorg.sigma_Rg
            guinier_autorg["Rg_error"].attrs["unit"] = radius_unit
            guinier_autorg["I0"] = autorg.I0
            guinier_autorg["I0_error"] = autorg.sigma_I0
            guinier_autorg["start_point"] = autorg.start_point
            guinier_autorg["end_point"] = autorg.end_point
            guinier_autorg["quality"] = autorg.quality
            guinier_autorg["aggregated"] = autorg.aggregated
            guinier_autorg["qₘᵢₙ·Rg"] =  autorg.Rg * res2.radial[autorg.start_point]
            guinier_autorg["qₘₐₓ·Rg"] =  autorg.Rg * res2.radial[autorg.end_point -1] 
            
    # Stage #4 Guinier plot generation:
        guinier = autorg or gpa #take one of the fit

        q, I, err = sasm.T[:3]
        mask = (I > 0) & numpy.isfinite(I) & (q > 0) & numpy.isfinite(q)
        if err is not None:
            mask &= (err > 0.0) & numpy.isfinite(err)
        mask = mask.astype(bool)
        if autorg:
            Rg = guinier.Rg
            I0 = guinier.I0
#             first_point = guinier.start_point
#             last_point = guinier.end_point
            intercept = numpy.log(I0)
            slope = -Rg * Rg / 3.0
            end = numpy.where(q > 1.5 / Rg)[0][0]
            mask[end:] = False

        q2 = q[mask] ** 2
        logI = numpy.log(I[mask])
        dlogI = err[mask] / logI
        q2_ds = guinier_data.create_dataset("q2", data=q2.astype(numpy.float32))
        q2_ds.attrs["unit"] = radius_unit+"⁻²"
        q2_ds.attrs["long_name"] = "q² (%s⁻²)"%radius_unit
        q2_ds.attrs["interpretation"] = "spectrum"
        lnI_ds = guinier_data.create_dataset("logI", data=logI.astype(numpy.float32))
        lnI_ds.attrs["long_name"] = "log(I)"
        lnI_ds.attrs["interpretation"] = "spectrum"
        erI_ds = guinier_data.create_dataset("errors", data=dlogI.astype(numpy.float32))
        erI_ds.attrs["interpretation"] = "spectrum"
                                             
        if autorg:
            guinier_data["fit"] = intercept + slope * q2
            guinier_data["fit"].attrs["slope"] = slope
            guinier_data["fit"].attrs["intercept"] = intercept
        
        guinier_data_attrs = guinier_data.attrs
        guinier_data_attrs["signal"] = "logI"
        guinier_data_attrs["axes"] = "q2"
        guinier_data_attrs["auxiliary_signals"] = "fit"
        guinier_grp.attrs["default"] = guinier_data.name
        if guinier is None:
            entry_grp.attrs["default"] = ai2_data.name
            self.log_error("No Guinier region found, data of dubious quality", do_raise=True)
    # Process 5: Kratky plot
        kratky_grp = nxs.new_class(entry_grp, "5_dimensionless_Kratky_plot", "NXprocess")
        kratky_grp["sequence_index"] = 5
        kratky_grp["program"] = "freesas.autorg"
        kratky_grp["version"] = freesas.version
        kratky_grp["date"] = get_isotime()
        kratky_data = nxs.new_class(kratky_grp, "results", "NXdata")
        kratky_grp.attrs["default"] = kratky_data.name
    # Stage #5 Kratky plot generation:
        Rg = guinier.Rg
        I0 = guinier.I0
        xdata = q * Rg
        ydata = xdata * xdata * I / I0
        dy = xdata * xdata * err / I0
        qRg_ds = kratky_data.create_dataset("qRg", data=xdata.astype(numpy.float32))
        qRg_ds.attrs["interpretation"] = "spectrum"
        qRg_ds.attrs["long_name"] = "q·Rg (unit-less)"
        k_ds = kratky_data.create_dataset("q2Rg2I/I0", data = ydata.astype(numpy.float32))
        k_ds.attrs["interpretation"] = "spectrum"
        k_ds.attrs["long_name"] = "q²Rg²I(q)/I₀"
        ke_ds = kratky_data.create_dataset("errors", data=dy.astype(numpy.float32))
        ke_ds.attrs["interpretation"] = "spectrum"
        kratky_data_attrs = kratky_data.attrs
        kratky_data_attrs["signal"] = "q2Rg2I/I0"
        kratky_data_attrs["axes"] = "qRg"
        
    # stage 6: Pair distribution function, what is the equivalent of datgnom
        bift_grp = nxs.new_class(entry_grp, "6_indirect_Fourier_transformation", "NXprocess")
        bift_grp["sequence_index"] = 6
        bift_grp["program"] = "freesas.bift"
        bift_grp["version"] = freesas.version
        bift_grp["date"] = get_isotime()
        bift_data = nxs.new_class(bift_grp, "results", "NXdata")
        cfg_grp = nxs.new_class(bift_grp, "configuration", "NXcollection")
    # Process stage6, i.e. perform the IFT
        try:
            bo = BIFT(q, I, err) 
            cfg_grp["Rg"] = guinier.Rg
            cfg_grp["npt"] = npt = 64 
            cfg_grp["Dmax÷Rg"] = 3
            Dmax = bo.set_Guinier(guinier, Dmax_over_Rg=3)
            # Pretty limited quality as we have real time constrains
            
            # First scan on alpha:
            cfg_grp["alpha_sup"] = alpha_max = bo.guess_alpha_max(npt)
            cfg_grp["alpha_inf"] = 1/alpha_max
            cfg_grp["alpha_scan_steps"] = 11
            
            key = bo.grid_scan(Dmax, Dmax, 1,
                               1.0 / alpha_max, alpha_max, 11, npt)
            Dmax, alpha = key[:2]
            # Then scan on Dmax:
            cfg_grp["Dmax_sup"] = guinier.Rg*4
            cfg_grp["Dmax_inf"] = guinier.Rg*2
            cfg_grp["Dmax_scan_steps"] = 5
            key = bo.grid_scan(guinier.Rg*2,guinier.Rg*4, 5,
                               alpha, alpha, 1, npt)
            Dmax, alpha = key[:2]
            if bo.evidence_cache[key].converged:
                bo.update_wisdom()
                use_wisdom = True
            else:
                use_wisdom = False
            res = minimize(bo.opti_evidence, (Dmax, log(alpha)), args=(npt, use_wisdom), method="powell")
            cfg_grp["Powell_steps"] = res.nfev
            cfg_grp["Monte-Carlo_steps"] = 0
        except Exception as error:
            bift_grp["Failed"] = "%s: %s"%(error.__class__.__name__, error)
            bo = None
        else:
            stats = bo.calc_stats()
            bift_grp["alpha"] = stats.alpha_avg
            bift_grp["alpha_error"] = stats.alpha_std
            bift_grp["Dmax"]=stats.Dmax_avg
            bift_grp["Dmax_error"]=stats.Dmax_std
            bift_grp["S0"]=stats.regularization_avg
            bift_grp["S0_error"]=stats.regularization_std
            bift_grp["Chi2r"]=stats.chi2r_avg
            bift_grp["Chi2r_error"]=stats.chi2r_std
            bift_grp["logP"]=stats.evidence_avg
            bift_grp["logP_error"]=stats.evidence_std
            bift_grp["Rg"]=stats.Rg_avg
            bift_grp["Rg_error"]=stats.Rg_std
            bift_grp["I0"]=stats.I0_avg
            bift_grp["I0_error"]=stats.I0_std
            #Now the plot:
            r_ds = bift_data.create_dataset("r", data=stats.radius.astype(numpy.float32))
            r_ds.attrs["interpretation"] = "spectrum"
            
            r_ds.attrs["unit"] = radius_unit
            r_ds.attrs["long_name"] = "radius r(%s)"%radius_unit
            p_ds = bift_data.create_dataset("p(r)", data=stats.density_avg.astype(numpy.float32))
            p_ds.attrs["interpretation"] = "spectrum"
            bift_data["errors"] = stats.density_std
            bift_data.attrs["signal"] = "p(r)"
            bift_data.attrs["axes"] = "r"
            
            r = stats.radius
            T = numpy.outer(q, r / pi)
            T = (4 * pi * (r[-1] - r[0]) / (len(r) - 1)) * numpy.sinc(T)
            bift_ds = ai2_data.create_dataset("BIFT", data=T.dot(stats.density_avg).astype(numpy.float32))
            bift_ds.attrs["interpretation"] = "spectrum"
            ai2_data.attrs["auxiliary_signals"] = "BIFT"
            bift_grp.attrs["default"] = bift_data.name
        #Finally declare the default entry and default dataset ...
        #overlay the BIFT fitted data on top of the scattering curve 
        entry_grp.attrs["default"] = ai2_data.name


    @staticmethod
    def read_nexus(filename):
        "return some NexusJuice from a HDF5 file "
        with Nexus(filename, "r") as nxsr:
            entry_grp = nxsr.get_entries()[0]
            h5path = entry_grp.name
            nxdata_grp = nxsr.h5[entry_grp.attrs["default"]]
            signal = nxdata_grp.attrs["signal"]
            axis = nxdata_grp.attrs["axes"]
            I = nxdata_grp[signal][()]
            q = nxdata_grp[axis][()]
            npt = len(q)
            unit = pyFAI.units.to_unit(axis+"_"+nxdata_grp[axis].attrs["units"])
            integration_grp = nxdata_grp.parent
            poni = str(integration_grp["configuration/file_name"][()]).strip()
            if not os.path.exists(poni):
                poni = str(integration_grp["configuration/data"][()]).strip()
            polarization = integration_grp["configuration/polarization_factor"][()]
            method = IntegrationMethod.select_method(**json.loads(integration_grp["configuration/integration_method"][()]))[0]
            instrument_grp = nxsr.get_class(entry_grp, class_type="NXinstrument")[0]
            detector_grp = nxsr.get_class(instrument_grp, class_type="NXdetector")[0]
            mask = detector_grp["pixel_mask"].attrs["filename"]
            mono_grp = nxsr.get_class(instrument_grp, class_type="NXmonochromator")[0]
            energy = mono_grp["energy"][()]
            img_grp = nxsr.get_class(entry_grp["3_time_average"], class_type="NXdata")[0]
            image2d = img_grp["intensity_normed"][()]
            error2d = img_grp["intensity_std"][()]
            sample_grp = nxsr.get_class(entry_grp, class_type="NXsample")[0]
            buffer = sample_grp["buffer"]
            concentration = sample_grp["concentration"]
        return NexusJuice(filename, h5path, npt, unit, q, I, poni, mask, energy, polarization, method, image2d, error2d, buffer, concentration)
