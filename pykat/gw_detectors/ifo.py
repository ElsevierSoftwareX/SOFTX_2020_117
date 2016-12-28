from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import math
import copy
import warnings
import cmath
from pykat import finesse
import pykat.components
import pykat.exceptions as pkex
import pykat.external.peakdetect as peak
import matplotlib.pyplot as plt
import pkg_resources
from scipy.optimize import fmin

global nsilica, clight
nsilica = 1.44963098985906
clight = 299792458.0

class aLIGO(object):
    """
    Object storing the Finesse input file and some auxilliary information
    about Avanced LIGO interferometers.

    References:
    [1] C. Bond `How to stay in shape: Overcoming beam and mirror distortions
        in advanced gravitational wave interferometers', PhD thesis, University
        of Birmingham, 2014, http://etheses.bham.ac.uk/5223/2/Bond14PhD.pdf
    [2] D. Martynov, `Lock Acquisition and Sensitivity Analysis of Advanced
        LIGO Interferometers', PhD thesis, Caltech, 2015
        http://thesis.library.caltech.edu/8899/1/DenisMartynovThesis.pdf
    [3] A. Staley, `Locking the Advanced LIGO Gravitational Wave Detector:
        with a focus on the Arm Length Stabilization Technique', PhD thesis,
        Columbia University, 2015
        https://academiccommons.columbia.edu/catalog/ac:189457
        
    """

    def __init__(self, _name="default", katfile=None, debug=False):
        names = ['default', 'LLO', 'LHO']
        if debug:
            self.kat = finesse.kat(tempdir=".",tempname="test")
        else:
            self.kat = finesse.kat()

        if katfile:
            self.kat.loadKatFile(katfile)
        else:
            if _name not in names: # TODO different files not yet implemented
                printf("aLIGO name `{}' not recognised, must be 'default', 'LLO' or 'LHO'",_name)
            _data_path=pkg_resources.resource_filename('pykat.gw_detectors','finesse_files/')
            #print(data_path)
            self.kat.loadKatFile(_data_path+"aLIGO.kat")

        # ----------------------------------------------------------------------
        # set variables to zero first
        self.DCoffset = 0.0
        self.DCoffsetW = 0.0

        # ----------------------------------------------------------------------
        # get and derive parameters from the kat file

        # get main frequencies
        if "f1" in self.kat.constants.keys():
            self.f1 = float(self.kat.constants["f1"].value)
        else:
            self.f1 = 9099471.0
        if "f2" in self.kat.constants.keys():
            self.f2 = float(self.kat.constants["f2"].value)
        else:
            self.f2 = 5.0 * self.f1
        if "f3" in self.kat.constants.keys():
            self.f3 = float(self.kat.constants["f3"].value)
        # TODO add else here!
        
        # defining a dicotionary for the main mirror positions (tunings)
        self.tunings = {}
        self.tunings = self.get_tunings(self.kat)
        # compute lengths such as PRC lentgth from individual lengths 
        self.compute_derived_lengths(self.kat)
        # check modultion frequencies
        if (5 * self.f1 != self.f2):
            print(" ** Warning: modulation frequencies do not match: 5*f1!=f2")
            
        # ----------------------------------------------------------------------
        # define ports and signals 
        
        # useful ports
        self.POP_f1  = port("POP_f1",   "nPOP",  self.f1, phase=101)
        self.POP_f2  = port("POP_f2",   "nPOP",  self.f2, phase=13)
        self.REFL_f1 = port("REFL_f1",  "nREFL", self.f1, phase=101)
        self.REFL_f2 = port("REFL_f2",  "nREFL", self.f2, phase=14)
        self.AS_DC   = port("AS_DC", "nSRM2")
        self.POW_BS  = port("PowBS", "nPRBS*")
        self.POW_X   = port("PowX",  "nITMX2")
        self.POW_Y   = port("PowY",  "nITMY2")

        # pretune DOF
        self.preARMX =  DOF("ARMX", self.POW_X,   "", "ETMX", 1, 1.0)
        self.preARMY =  DOF("ARMY", self.POW_Y,   "", "ETMY", 1, 1.0)
        self.preMICH =  DOF("AS"  , self.AS_DC,   "", ["ITMX", "ETMX", "ITMY", "ETMY"], [1,1,-1,-1], 6.0)
        self.prePRCL =  DOF("PRCL", self.POW_BS,  "", "PRM",  1, 10.0)
        self.preSRCL =  DOF("SRCL", self.AS_DC,   "", "SRM",  1, 10.0)
        
        # control scheme as in [1] Table C.1  
        self.PRCL =  DOF("PRCL", self.POP_f1,  "I", "PRM", 1, 100.0)
        self.MICH =  DOF("MICH", self.POP_f2,  "Q", ["ITMX", "ETMX", "ITMY", "ETMY"], [1,1,-1,-1], 100.0)
        self.CARM =  DOF("CARM", self.REFL_f1, "I", ["ETMX", "ETMY"], [1, 1], 1.5)
        self.DARM =  DOF("DARM", self.AS_DC,   "",  ["ETMX", "ETMY"], [1,-1], 1.0)
        self.SRCL =  DOF("SRCL", self.REFL_f2, "I", "SRM", 1, 1e2)


            
    def adjust_PRC_length(self, kat):
        """
        Adjust PRC length so that it fulfils the requirement
        lPRC = (N+1/2) * c/(2*f1), see [1] equation C.1
        In the current design N=3
        """
        print("-- adjusting PRC length")
        ltmp = 0.5 * clight / self.f1
        delta_l = 3.5 * ltmp - self.lPRC
        print("  adusting kat.lp1.L by {}m".format(delta_l))
        kat.lp1.L += delta_l
        self.compute_derived_lengths(kat)

    def check_f1_PRC_resonance(self, _kat):
        """
        Plot the sideband amplitudes for mdulation frequecy
        f1 (~ 9MHz) in the PRC, to check the resonance
        condition.
        """
        kat = _kat.deepcopy()
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        startf = self.f1-400.0
        stopf  = self.f1+400.0
        code = """
        ad f1p {0} nPRM2
        ad f1m -{0} nPRM2
        xaxis mod1 f lin {1} {2} 200
        put f1p f $x1 
        put f1m f $mx1 
        """.format(self.f1, startf, stopf)
        kat.parseCommands(code)
        out = kat.run()
        ax.plot(out.x-self.f1,np.abs(out["f1p"]), label=" f1")
        ax.plot(out.x-self.f1,np.abs(out["f1m"]), label="-f1")
        ax.set_xlim([np.min(out.x-self.f1), np.max(out.x-self.f1)])
        ax.set_xlabel("delta_f1 [Hz]")
        ax.set_ylabel('sqrt(W) ')
        ax.grid()
        ax.legend()
        plt.tight_layout()
        plt.show(block=0)

        
    def compute_derived_lengths(self, kat, verbose=False):
        """
        Compute derived length from individual space components.
        Design values are currently:
        lPRC = 57.656, lSRC = 56.008, lSchnupp = 0.08
        and the individual lengths:
        PRC: Lp1 16.6107, Lp2 16.1647, Lp3 19.5381 
        SRC: Ls1 15.7586, Ls2 15.4435, Ls3 19.3661
        """
        # distances between HR surfaces:
        self.lpr = kat.lp1.L + kat.lp2.L + kat.lp3.L
        self.lx = kat.lx1.L + kat.BSsub1.L * kat.BSsub1.n + kat.ITMXsub.L * kat.ITMXsub.n
        self.ly = kat.ly1.L + kat.ITMYsub.L * kat.ITMYsub.n
        self.lsr = kat.ls1.L + kat.ls2.L + kat.ls3.L + kat.BSsub2.L * kat.BSsub2.n
        # resulting combined distances (single, not roundtrip)
        self.lMI =  0.5 * (self.lx + self.ly)
        self.lPRC = self.lpr + self.lMI
        self.lSRC = self.lsr + self.lMI
        self.lSchnupp = self.lx - self.ly
        if verbose:
            print("-- small MI lengths")
            print(" lx = {}m, ly = {}m".format(np.round(self.lx, 4), np.round(self.ly, 4)))
            #print(" lpr = {}m, lsr = {}m".format(np.round(self.lpr, 4),  np.round(self.lsr),4))
            #print(" lMI = {}m, lSchnupp = {}m".format(np.round(self.lMI, 4) , self.lSchnupp))
            print(" lSchnupp = {}m".format(np.round(self.lSchnupp,4)))
            print(" lPRC = {}m, lSRC= {}m".format(self.lPRC, self.lSRC))
        
    def get_tunings(self, kat):
        self.tunings["maxtem"] = kat.maxtem
        self.tunings["PRM"]    = kat.PRM.phi
        self.tunings["ITMX"]   = kat.ITMX.phi
        self.tunings["ETMX"]   = kat.ETMX.phi
        self.tunings["ITMY"]   = kat.ITMY.phi
        self.tunings["ETMY"]   = kat.ETMY.phi
        self.tunings["BS"]     = kat.BS.phi
        self.tunings["SRM"]    = kat.SRM.phi
        return self.tunings
    
    def set_tunings(self, kat, tunings):
        kat.maxtem   = tunings["maxtem"]
        kat.PRM.phi  = tunings["PRM"]  
        kat.ITMX.phi = tunings["ITMX"]  
        kat.ETMX.phi = tunings["ETMX"]  
        kat.ITMY.phi = tunings["ITMY"]  
        kat.ETMY.phi = tunings["ETMY"]  
        kat.BS.phi   = tunings["BS"]  
        kat.SRM.phi  = tunings["SRM"]  

    def apply_lock_feedback(self, kat, out):
        tuning = self.get_tunings(kat)
        if "ETMX_lock" in out.keys():
            tuning["ETMX"] += float(out["ETMX_lock"])
        else:
            print(" ** Warning: could not find ETMX lock")
        if "ETMY_lock" in out.keys():
            tuning["ETMY"] += float(out["ETMY_lock"])
        else:
            print(" ** Warning: could not find ETMY lock")
        if "PRCL_lock" in out.keys():
            tuning["PRM"]  += float(out["PRCL_lock"])
        else:
            print(" ** Warning: could not find PRCL lock")
        if ("MICH_lock" in out.keys()) and ("ITMY_lock" in out.keys()):
            tuning["ITMX"] += float(out["MICH_lock"])
            tuning["ITMY"] += float(out["ITMY_lock"])
        else:
            print(" ** Warning: could not find MICH (ITMY) lock")
        if "SRCL_lock" in out.keys():
            tuning["SRM"]  += float(out["SRCL_lock"])
        else:
            print(" ** Warning: could not find SRCL lock")
        self.set_tunings(kat, tuning)

        
    def pretune(self, _kat, pretune_precision=1.0e-4):
        print("-- pretuning interferometer to precision {0:2g} deg = {1:2g} m".format(pretune_precision, pretune_precision*_kat.lambda0/360.0))
        kat=_kat.deepcopy()
        #_kat.ITMX.phi=0.0
        #_kat.ITMY.phi=0.0
        print("  scanning X arm (maximising power)")
        kat1 = kat.deepcopy()
        make_transparent(kat1,["PRM","SRM"])
        make_transparent(kat1,["ITMY", "ETMY"])
        kat1.BS.setRTL(0.0,1.0,0.0) # set BS refl. for X arm
        phi, precision = self.scan_to_precision(kat1, self.preARMX, pretune_precision)
        phi=round(phi/pretune_precision)*pretune_precision
        phi=round_to_n(phi,5)
        print("  found max/min at: {} (precision = {:2g})".format(phi, precision))
        self.preARMX.apply_tuning(_kat,phi)
    
        print("  scanning Y arm (maximising power)")
        kat = _kat.deepcopy()
        make_transparent(kat,["PRM","SRM"])
        make_transparent(kat,["ITMX", "ETMX"])
        kat.BS.setRTL(1.0,0.0,0.0) # set BS refl. for Y arm
        phi, precision = self.scan_to_precision(kat, self.preARMY, pretune_precision)
        phi=round(phi/pretune_precision)*pretune_precision
        phi=round_to_n(phi,5)
        print("  found max/min at: {} (precision = {:2g})".format(phi, precision))
        self.preARMY.apply_tuning(_kat,phi)
    
        print("  scanning MICH (minimising power)")
        kat = _kat.deepcopy()
        make_transparent(kat,["PRM","SRM"])
        phi, precision = self.scan_to_precision(kat, self.preMICH, pretune_precision, minmax="min", precision=30.0)
        phi=round(phi/pretune_precision)*pretune_precision
        phi=round_to_n(phi,5)
        print("  found max/min at: {} (precision = {:2g})".format(phi, precision))
        self.preMICH.apply_tuning(_kat,phi, add=True)

        print("  scanning PRCL (maximising power)")
        kat = _kat.deepcopy()
        make_transparent(kat,["SRM"])
        phi, precision = self.scan_to_precision(kat, self.prePRCL, pretune_precision)
        phi=round(phi/pretune_precision)*pretune_precision
        phi=round_to_n(phi,5)
        print("  found max/min at: {} (precision = {:2g})".format(phi, precision))
        self.prePRCL.apply_tuning(_kat,phi)

        print("  scanning SRCL (maximising carrier power, then adding 90 deg)")
        kat = _kat.deepcopy()
        phi, precision = self.scan_to_precision(kat, self.preSRCL, pretune_precision, phi=0)
        phi=round(phi/pretune_precision)*pretune_precision
        phi=round_to_n(phi,4)-90.0
        print("  found max/min at: {} (precision = {:2g})".format(phi, precision))
        self.preSRCL.apply_tuning(_kat,phi)
        
    def scan_to_precision(self, kat, DOF, pretune_precision, minmax="max", phi=0.0, precision=60.0):
        while precision>pretune_precision*DOF.scale:
            out = scan_DOF(kat, DOF, xlimits = [phi-1.5*precision, phi+1.5*precision])
            phi, precision = find_peak(out, DOF.port.portName, minmax=minmax)
            #print("** phi= {}".format(phi))
        return phi, precision

    def pretune_status(self, _kat):
        kat = _kat.deepcopy()
        kat.verbose = False
        kat.noxaxis = True
        pretune_DOFs = [self.preARMX, self.preARMY, self.prePRCL, self.preMICH, self.preSRCL]
        _detStr=""
        for p in pretune_DOFs:
            _sigStr = p.port.signal(kat)
            _detStr = "\n".join([_detStr, _sigStr])
        kat.parseCommands(_detStr)
        out = kat.run()
        Pin = float(kat.L0.P)
    
        tunings = self.get_tunings(kat)
        print(" .------------------------------------------------------.")
        print(" | pretuned for maxtem = {:4}          (-1 = 'off')     |".format(tunings["maxtem"]))
        keys_t = list(tunings.keys())
        keys_t.remove("maxtem")
        print(" .------------------------------------------------------.")
        print(" | port   power[W] pow. ratio | optics   tunings        |")
        print(" +----------------------------|-------------------------+")
        idx_p = 0
        idx_t = 0
        run_p = True
        run_t = True
        while (run_p or run_t):
            if idx_p < len(pretune_DOFs):
                p = pretune_DOFs[idx_p]
                print(" | {:5}: {:8.3g} {:8.3g}   |".format(p.name, float(out[p.port.name]), float(out[p.port.name])/Pin),end="")
                idx_p +=1
            else:
                print(" |                            |", end="")
                run_p = False
            if idx_t < len(keys_t):
                t=keys_t[idx_t]
                print(" {:5}: {:9.3g}        |".format(t, float(self.tunings[t])))
                idx_t +=1
            else:
                print("                         |")
                run_t = False
        print(" `------------------------------------------------------'")
                     
    def power_ratios(self, _kat):
        kat = _kat.deepcopy()
        kat.verbose = False
        kat.noxaxis = True
    
        ports = [self.POW_X, self.POW_Y, self.AS_DC, self.POW_BS]
        _detStr = ""
        for p in ports:
            _sigStr = p.signal(kat)
            _detStr = "\n".join([_detStr, _sigStr])
        kat.parseCommands(_detStr)
        out = kat.run()
        
        Pin = float(kat.L0.P)
    
        print("-- power ratios (Pin = {0:.3g} W)".format(Pin))
        for p in ports:
            print(" {0:6} = {1:8.3g} W ({0:6}/Pin = {2:8.2g})" .format(p.name, float(out[p.name]), float(out[p.name])/Pin))
        return

    def plot_pretuning_powers(self, _kat, xlimits=[-10,10]):
        kat = _kat.deepcopy()
        kat.verbose = False
        kat.noxaxis = True
        dofs = [self.preARMX, self.preARMY, self.preMICH, self.prePRCL]
        idx=1
        fig = plt.figure()
        for d in dofs:
            ax = fig.add_subplot(2,2,idx)
            idx+=1
            out = scan_DOF(kat, d, xlimits = np.multiply(d.scale, xlimits), relative = True)
            ax.semilogy(out.x,out[d.signal_name(kat)])
            ax.set_xlim([np.min(out.x), np.max(out.x)])
            ax.set_xlabel("phi [deg] {}".format(d.optics[0]))
            ax.set_ylabel('{} [W] '.format(d.signal_name(kat)))
            ax.grid()
        plt.tight_layout()
        plt.show(block=0)

    def plot_error_signals(self, _kat, xlimits=[-1,1]):
        kat = _kat.deepcopy()
        kat.verbose = False
        kat.noxaxis = True
        dofs = [self.DARM, self.CARM, self.PRCL, self.SRCL,  self.MICH]
        idx=1
        fig = plt.figure()
        for d in dofs:
            ax = fig.add_subplot(2,3,idx)
            idx+=1
            out = scan_DOF(kat, d, xlimits = np.multiply(d.scale,xlimits), relative = True)
            ax.plot(out.x,out[d.signal_name(kat)])
            ax.set_xlim([np.min(out.x), np.max(out.x)])
            ax.set_xlabel("{} [deg]".format(d.name))
            ax.set_ylabel('{} [W] '.format(d.port.name))
            ax.grid()
        plt.tight_layout()
        plt.show(block=0)

    def find_DC_offset(self, _kat, AS_power, precision=1e-4):
        """
        Returns the DC offset of DARM that corrponds to the
        specified power in the AS power.
        """
        print("-- finding DC offset for AS power of {:3g} W".format(AS_power))
        kat = _kat.deepcopy()
        kat.verbose = False
        kat.noxaxis = True
        _sigStr = self.AS_DC.signal(kat)
        kat.parseCommands(_sigStr)
        Xphi = float(kat.ETMX.phi)
        Yphi = float(kat.ETMY.phi)

        def powerDiff(phi, kat, Xphi, Yphi, AS_power):
            kat.ETMY.phi = Yphi + phi
            kat.ETMX.phi = Xphi - phi
            out=kat.run()
            #print(out[self.AS_DC.name]-AS_power)
            return np.abs(out[self.AS_DC.name]-AS_power)

        out=fmin(powerDiff,0,xtol=precision,ftol=1e-3,args=(kat, Xphi, Yphi, AS_power))
        print("  DC offset for AS_DC={} W is: {}".format(AS_power, out[0]))
        self.DCoffset = round(out[0],6)
        self.DCoffsetW = AS_power
        return self.DCoffset

    def generate_errsig_block(self, kat, noplot=False):
        sigDARM = self.DARM.signal(kat)
        sigCARM = self.CARM.signal(kat)
        sigPRCL = self.PRCL.signal(kat)
        sigMICH = self.MICH.signal(kat)
        sigSRCL = self.SRCL.signal(kat)
        code1 = "\n".join([sigDARM, sigCARM, sigPRCL, sigMICH, sigSRCL])

        if noplot:
            nameDARM = self.DARM.signal_name(kat)
            nameCARM = self.CARM.signal_name(kat)
            namePRCL = self.PRCL.signal_name(kat)
            nameMICH = self.MICH.signal_name(kat)
            nameSRCL = self.SRCL.signal_name(kat)
            code2 = """
            noplot {}
            noplot {}
            noplot {}
            noplot {}
            noplot {}
            """.format(nameDARM, nameCARM, namePRCL, nameMICH, nameSRCL)
            return "".join([code1, code2])
        else:
            return code1
            
        
    def generate_lock_block(self, _kat, _gains=None, _accuracies=None, verbose=False):
        """
        gains: optical gain is in W per rad
        accuracies: error signal threshold in W
        rms, estimated loop noise rms m

        to compute accuracies from rms, we convert
        rms to radians as rms_rad = rms * 2 pi/lambda
        and then multiply by the optical gain.
        """
        kat = _kat.deepcopy()
        if _gains == None:
            ogDARM = optical_gain(kat, self.DARM, self.DARM)
            ogCARM = optical_gain(kat, self.CARM, self.CARM)
            ogPRCL = optical_gain(kat, self.PRCL, self.PRCL)
            ogMICH = optical_gain(kat, self.MICH, self.MICH)
            ogSRCL = optical_gain(kat, self.SRCL, self.SRCL)
            if verbose:
                print("-- optical gains:")
                print("  DARM: {}".format(ogDARM))
                print("  CARM: {}".format(ogCARM))
                print("  PRCL: {}".format(ogPRCL))
                print("  MICH: {}".format(ogMICH))
                print("  SRCL: {}".format(ogSRCL))
            gains = [ ogDARM, ogCARM, ogPRCL, ogMICH, ogSRCL]
        else:
            gains = _gains.copy()

        rms = [1e-13, 1e-12, 1e-11, 1e-11, 1e-11]
        factor = 2.0 * math.pi / kat.lambda0
        if _accuracies == None:
            accDARM = round_to_n(np.abs(factor * rms[0] * gains[0]),2) 
            accCARM = round_to_n(np.abs(factor * rms[1] * gains[1]),2) * 0.1 # manually tuned
            accPRCL = round_to_n(np.abs(factor * rms[2] * gains[2]),2) * 0.1 # manually tuned
            accMICH = round_to_n(np.abs(factor * rms[3] * gains[3]),2)
            accSRCL = round_to_n(np.abs(factor * rms[4] * gains[4]),2) * 50.0 # manually tuned
            acc = [accDARM, accCARM, accPRCL, accMICH, accSRCL]
        else:
            acc = _accuracies.copy()

        nameDARM = self.DARM.signal_name(kat)
        nameCARM = self.CARM.signal_name(kat)
        namePRCL = self.PRCL.signal_name(kat)
        nameMICH = self.MICH.signal_name(kat)
        nameSRCL = self.SRCL.signal_name(kat)
        code1 = """
        %%% FTblock locks
        ###########################################################################
        set AS_f2_I_re {} re
        set CARM_err {} re
        set PRCL_err {} re
        set MICH_err {} re
        set SRCL_err {} re
        func DARM_err = $AS_f2_I_re - {}
        """.format(nameDARM, nameCARM, namePRCL, nameMICH, nameSRCL, self.DCoffsetW)

        factor = 0.4 * -1.0 * 180 / math.pi # 0.2 because of multiple locks cross talk
        gainDARM = round_to_n(factor / ogDARM, 2)
        gainCARM = round_to_n(0.01 * factor / ogCARM, 2) # factor 0.01 for better gain hirarchy
        gainPRCL = round_to_n(2.0  * factor / ogPRCL, 2) # manually tuned
        gainMICH = round_to_n(1.0  * factor / ogMICH, 2) # manually tuned
        gainSRCL = round_to_n(0.05 * factor / ogSRCL, 2) # gain hirrchy with MICH
        
        code2 = """
        lock DARM_lock $DARM_err {} {}
        lock CARM_lock $CARM_err {} {} 
        lock PRCL_lock $PRCL_err {} {}
        lock MICH_lock $MICH_err {} {} 
        lock SRCL_lock $SRCL_err {} {} 
        """.format(gainDARM, acc[0], gainCARM, acc[1], gainPRCL, acc[2], gainMICH, acc[3], gainSRCL, acc[4])
        
        code3 = """
        noplot ITMY_lock
        func ITMY_lock = (-1.0) * $MICH_lock 
        func ETMX_lock = $CARM_lock + $MICH_lock + $DARM_lock
        func ETMY_lock = $CARM_lock - $MICH_lock - $DARM_lock

        put* PRM     phi     $PRCL_lock
        put* ITMX    phi     $MICH_lock
        put* ITMY    phi     $ITMY_lock
        put* ETMX    phi     $ETMX_lock
        put* ETMY    phi     $ETMY_lock
        put* SRM     phi     $SRCL_lock

        noplot PRCL_lock
        noplot SRCL_lock
        noplot MICH_lock
        noplot DARM_lock
        noplot CARM_lock
        noplot ETMX_lock
        noplot ETMY_lock

        ###########################################################################
        %%% FTend locks
        """
        return "".join([code1, code2, code3])

# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
class DOF(object):
    """
    Defining a degree of freedom for the interferometer, includes the
    objects and how to move them, and the default output port to read
    out the DOF signal.
    """
    def __init__(self, _DOFName, _port, _quad, _optics, _factors, _scale, sigtype="z"):
        self.name = _DOFName
        self.port = _port
        self.quad = _quad
        self.sigtype = sigtype
        self.optics=make_list_copy(_optics)
        self.factors=make_list_copy(_factors)
        # scaling factor, to compensate for lower sensitivity compared
        # to DARM (in tuning plots for example)
        # Thus DARM has a scale of 1, all other DOFs a scale >1
        self.scale = _scale

    def apply_tuning(self, kat, phi, add=False):
        for idx, o in enumerate(self.optics):
            if add:
                kat.components[o].phi += phi * self.factors[idx]
            else:
                kat.components[o].phi = phi * self.factors[idx]
            
    def signal(self, kat):
        return self.port.signal(kat, self.quad, sigtype=self.sigtype)
    def signal_name(self, kat):
        return self.port.signal_name(kat, self.quad, sigtype=self.sigtype)

    def transfer(self, kat, fsig, phase2=None):
        return self.port.transfer(kat, self.quad, fsig=fsig, phase2=phase2, sigtype=self.sigtype)
    def transfer_name(self, kat):
        return self.port.transfer_name(kat, self.quad)

    def fsig(self, _fsigName, fsig=1.0):
        _fsigStr= ""
        for idx, o in enumerate(self.optics):
            phase = 0.0
            if self.factors[idx] == -1:
                phase = 180.0
            _fsigStr = "\n".join([_fsigStr, "fsig {} {} {} {} ".format(_fsigName, o, fsig, phase)])
        return _fsigStr

# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
class port(object):
    """
    Defining an output port for the interferometer, can be either a
    pd or a pd1 detector (for error signal generation).
    """
    def __init__(self, _portName, _nodeNames, f=None, phase=None):
        self.portName = _portName
        self.nodeNames = make_list_copy(_nodeNames)
        self.f=f            # demodulation frequency, float
        self.phase = phase  # demodulation frequency for I quadrature, float
        self.name = self.portName            
    
    def check_nodeName(self, kat):
        self.nodeName = None
        for node in self.nodeNames:
            _node = node
            if _node[-1] == "*":
                _node = _node[:-1]
            if _node in kat.nodes:
                self.nodeName=node
                break
        if self.nodeName==None:
            raise pkex.BasePyKatException("port {}: cannot find any of these nodes: '{}'".format(self.name,self.nodeNames))

    def amplitude_name(self, kat, f, n=None, m=None, sigtype="z"):
        name = self.name + "_ad"
        return name
    
    def amplitude(self, kat, f, n=None, m=None, sigtype="z"):
        self.check_nodeName(kat)
        name = self.amplitude_name(kat, f, n=n, m=m, sigtype=sigtype)
        if n==None and m==None:
            return "ad {} {} {}".format(name, f, self.nodeName)
        else:
            return "ad {} {} {} {} {}".format(name, f, n, m, self.nodeName)
    
    def signal_name(self, kat, quad="I", sigtype="z"):
        name = self.name
        if self.f!=None:
            name = self.name+"_"+quad
        return name
    
    def signal(self, kat, quad="I", sigtype="z"):
        self.check_nodeName(kat)
        name = self.signal_name(kat, quad=quad, sigtype=sigtype)
        if sigtype != "z":
                raise pkex.BasePyKatException("alignment signals are not implemented yet")            
        if self.f==None:
            return "pd {} {}".format(name, self.nodeName)
        else:
            if quad !="I" and quad != "Q":
                raise pkex.BasePyKatException("quadrature must be 'I' or 'Q'")            
            phase = self.IQ_phase(quad, self.phase)
            return "pd1 {} {} {} {}".format(name, self.f, phase, self.nodeName)
        
    def IQ_phase(self, quad, phase):
        if quad== "Q":
            phase = phase + 90.0
            if phase >=360.0 :
                phase -= 360.0
        return phase
        
    def transfer_name(self, kat, quad="I"):
        name = self.name
        if self.f!=None:
            name = self.name+"_"+quad
        return name
        
    def transfer(self, kat, quad="I", fsig=1.0, phase2=None, sigtype="z"):
        self.check_nodeName(kat)
        name = self.transfer_name(kat, quad=quad)
        if sigtype!="z":
                raise pkex.BasePyKatException("alignment signals are not implemented yet")            
        if self.f==None:
            if phase2 == None:
                return "pd1 {} {} {}".format(name, fsig, self.nodeName)
            else:
                return "pd1 {} {} {} {}".format(name, fsig, phase2, self.nodeName)
        else:
            if quad !="I" and quad != "Q":
                raise pkex.BasePyKatException("quadrature must be 'I' or 'Q'")            
            phase = self.IQ_phase(quad, self.phase)
            if phase2 == None:
                return "pd2 {} {} {} {} {}".format(name, self.f , phase, fsig, self.nodeName)
            else:
                return "pd2 {} {} {} {} {} {}".format(name, self.f, phase, fsig, phase2, self.nodeName)
                                
        
def scan_optics_string(_optics, _factors, _varName, linlog="lin", xlimits=[-100, 100], steps=200, axis=1,relative=False):
    optics=make_list_copy(_optics)
    factors=make_list_copy(_factors)
    if len(optics) != len(factors):
        raise pkex.BasePyKatException("you must provide a factor for each optics")

    if linlog not in ["lin", "log"]: 
        raise pkex.BasePyKatException("linlog must be 'lin' or 'log'")
    _tuneStr  = "var {} 0\n".format(_varName)
    if axis==1:
        _tuneStr += "xaxis {} phi {} {} {} {}".format(_varName, linlog, xlimits[0], xlimits[1], steps)
    elif (axis==2 or axis==3): 
        _tuneStr += "x{}axis {} phi {} {} {} {}".format(axis, _varName, linlog, xlimits[0], xlimits[1], steps)
    else:
        raise pkex.BasePyKatException("axis must be 1, 2 or 3")
    _putStr = ""
    for idx, o in enumerate(optics):
        if factors[idx] == 1:
            _xStr="$x"
        elif factors[idx] == -1:
            _xStr="$mx"
        else:
            raise pkex.BasePyKatException("optics factors must be 1 or -1")
        if (relative):
            _putCmd = "put*"
        else:
            _putCmd = "put"
                
        _putStr = "\n".join([_putStr, "{} {} phi {}{}".format(_putCmd, o, _xStr, axis)])            
        _tuneStr += _putStr
    return _tuneStr
        
def make_transparent(kat, _components):
    """
    Function to make certain mirror or beamsplitter objects transparent
    Usage:
    make_transparent(kat, "PRM")
    or
    make_transparent(kat, ["PRM", "SRM"])
    The function modifies the kat object. Use something like 
    kat = copy.deepcopy(basekat)
    before if you want to preseve the original kat object
    """
    components=make_list_copy(_components)
    for o in kat.components.values():
        if o.name in components:
            if isinstance(o, pykat.components.AbstractMirrorComponent): 
                o.setRTL(0.0,1.0,0.0)
                components = [c for c in components if c != o.name]
            else:
                raise pkex.BasePyKatException("{} is not a mirror or beamsplitter".format(o))
            #print('  made {} transparent'.format(o))
    if len(components) != 0:
        raise pkex.BasePyKatException("Cannot find component {}".format(components))
    return kat

def reconnect_nodes(kat, component1, idx1, node_name):
    c_string = component1.getFinesseText()
    c = c_string[0].split()
    new_string = " ".join(c[:-2])
    nodes = {}
    nodes[0] = c[-2]
    nodes[1] = c[-1]
    nodes[idx1]=node_name
    new_string = new_string + " " + nodes[0] + " " + nodes[1]
    #print(" new string ='{}'".format(new_string))
    kat.parseCommands(new_string)

def remove_commands(kat, _commands):
    commands=make_list_copy(_commands)
    # removing commands
    for o in kat.commands.values():
        if o.name in commands:
            o.remove()
            commands = [c for c in commands if c != o.name]
            #print('  {} removed'.format(o))
    if len(commands) != 0:
        raise pkex.BasePyKatException("Cannot find command(s) {}".format(commands))
    return kat
    
def remove_components(kat, _components, component_in=None, component_out=None):
    components=make_list_copy(_components)
    if  kat.components[components[-1]].nodes[1]:
        node_in  = kat.components[components[-1]].nodes[1].name
    else:
        node_in = None
    node_out = kat.components[components[0]].nodes[0].name
    # removing components
    for o in kat.components.values():
        if o.name in components:
            o.remove()
            components = [c for c in components if c != o.name]
            #print('  {} removed'.format(o))
    if len(components) != 0:
        raise pkex.BasePyKatException("Cannot find component(s) {}".format(components))
    # reconnecting nodes if requested
    if component_in:
        reconnect_nodes(kat, kat.components[component_in],1, node_in)
    if component_out:
        reconnect_nodes(kat, kat.components[component_out],0, node_out)
    return kat
        
def BS_optical_path(thickness, n=nsilica, angle=45.0):
    """
    Compute optical path length in BS substrate, default
    parameters assume angle of incidence of 45 deg and fused
    silica substrate.
    thickness: substrate thickness [m]
    n: index of refraction of substrate
    angle: angle of incidence (in vacuum) [deg]
    """
    
    angle_subst = math.asin(math.sin(math.radians(angle))/n)
    L = thickness / math.cos(angle_subst) 
    
    return math.degrees(angle_subst), L

def scan_DOF(_kat, DOF, xlimits=[-100, 100], steps=200, relative=False): 
    kat = _kat.deepcopy()
    scan_string = scan_optics_string(DOF.optics, DOF.factors, "scan", linlog="lin", xlimits=xlimits, steps=steps, axis=1,relative=relative)
    kat.parseCommands(scan_string)
    sigStr = DOF.signal(kat)
    kat.parseCommands(sigStr)
    out = kat.run()
    return out

def scan_optics(_kat, _optics, _factors, xlimits=[-100, 100], steps=200,relative=False): 
    """
    Scans one or more optics (by changing its tuning).
    Parameters:
    optics: list of names of components to be tuned
    factors: list of scaling factors for the tuning for each optics, first element must be 1.0
    xlimits: limits of scan, defaults to [-100, 100]
    steps: number of steps to use in scan, default is 200
    Usage:
    scan_optis(kat, "PRM", 1)
    scan_optis(kat, ["ETMX", "ETMY", [1, -1])
    """
    kat = _kat.deepcopy()
    optics=make_list_copy(_optics)
    factors=make_list_copy(_factors)
    scan_string = scan_optics_string(optics, factors, "scan", linlog="lin", xlimits=xlimits, steps=steps, axis=1,relative=relative)
    kat.parseCommands(scan_string)
    out = kat.run()
    return out

def optical_gain(_kat, DOF_sig, DOF_det, f=10.0):
    kat = _kat.deepcopy()
    _fsigStr = DOF_sig.fsig("sig1", fsig=f)
    _detStr = DOF_det.transfer(kat, fsig=f)
    _detName = DOF_det.transfer_name(kat)
    kat.parseCommands(_fsigStr)
    kat.parseCommands(_detStr)
    kat.noxaxis = True
    kat.parseCommands("yaxis lin abs:deg")
    out = kat.run()
    return np.real(out[_detName])

def find_peak(out, detector, minmax='max', debug=False): 
    """
    Expects an output of a kat scan and find the max/min output on
    a sepcific detector. Useful for pre-tuning cavities or interferometers.
    Returns the tuning of the maximum (or minimum) and the precision
    Parameters:
    detector: name of detetor to use
    minmax, string to indicate maximum or minum detection, default 'max'
    Usage:
    find_peak(out, "pdout")
    find_peak(out, "pdout", minmax='min')
    """
    stepsize = out.x[1]-out.x[0]
    print("  stepsize (precision) of scan: {0:g}".format(stepsize))

    _max, _min = peak.peakdetect( out[detector],out.x, 1)
    
    if debug==True:
        plt.figure()
        plt.plot(out.x,out[detector])
        print("max: ")
        print(_max)
        print("min: ")
        print(_min)
        
    if minmax == 'max':
        X = [p[0] for p in _max]
        Y = [p[1] for p in _max]
        X_out = X[np.argmax(Y)]
        Y_out = np.max(Y)
    elif minmax == 'min':
        X = [p[0] for p in _min]
        Y = [p[1] for p in _min]
        X_out = X[np.argmin(Y)]
        Y_out = np.min(Y)
    else:
        raise pkex.BasePyKatException("maxmin must be 'max' or 'min'")
        
    if debug==True:
        plt.plot(X_out,Y_out,'o')
        plt.xlabel('tuning [deg]')
        plt.ylabel('{0} output'.format(detector))
        plt.show()
    return X_out, stepsize

def make_list_copy(_l):
    """
    Utility function, takes a list of strings or single string
    and returns a copy of a list, e.g.
    "string" copy to ["string"]
    ["string1", "string2"] copy to ["string1", "string2"]
    """
    if not isinstance(_l, list):
        _l = [_l]
    return _l[:] # copy the list, just to be save

def round_to_n(x, n):
    if not x: return 0
    power = -int(math.floor(math.log10(abs(x)))) + (n - 1)
    factor = (10 ** power)
    return round(x * factor) / factor
