from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import math
import copy
import warnings
import cmath
import inspect
import six 

from pykat import finesse
from pykat.finesse import BlockedKatFile

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

    def __init__(self, _name="default", katfile=None, verbose = False, debug=False):
        names = ['default', 'LLO', 'LHO']
        if debug:
            self.kat = finesse.kat(tempdir=".",tempname="test")
        else:
            self.kat = finesse.kat()
        self.kat.verbose=verbose
        self._data_path=pkg_resources.resource_filename('pykat.gw_detectors','finesse_files/')

        self.rawBlocks = BlockedKatFile()
        if katfile:
            self.kat.loadKatFile(katfile)
            self.rawBlocks.read(katfile)
        else:
            """
            if _name not in names: # TODO different files not yet implemented
                printf("aLIGO name `{}' not recognised, must be 'default', 'LLO' or 'LHO'",_name)
            """
            if _name != "default":
                printf("aLIGO name `{}' not recognised, using 'default'",_name)                
            self.kat.loadKatFile(self._data_path+"aLIGO.kat")
            self.rawBlocks.read(self._data_path+"aLIGO.kat")

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
        # check modultion frequencies
        if (5 * self.f1 != self.f2):
            print(" ** Warning: modulation frequencies do not match: 5*f1!=f2")
        
        # defining a dicotionary for the main mirror positions (tunings),
        # keys should include maxtem, phase and all main optics names
        self.tunings = dict.fromkeys(["maxtem", "phase", "PRM", "ITMX", "ETMX", "ITMY", "ETMY", "BS", "SRM"])
        self.tunings = get_tunings(self.kat, self.tunings)
        self.compute_derived_lengths(self.kat)
            
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
        
        self.__DOFs = {}
        
        for _ in inspect.getmembers(self, lambda x: isinstance(x, DOF)):
            self.__DOFs[_[0]] = _[1]

        self.lockNames = None

    @property
    def DOFs(self):
        return copy.copy(self.__DOFs)
    
    def adjust_PRC_length(self, kat, verbose=False):
        """
        Adjust PRC length so that it fulfils the requirement
        lPRC = (N+1/2) * c/(2*f1), see [1] equation C.1
        In the current design N=3
        """
        vprint(verbose, "-- adjusting PRC length")
        ltmp = 0.5 * clight / self.f1
        delta_l = 3.5 * ltmp - self.lPRC
        vprint(verbose, "   adusting kat.lp1.L by {:.4g}m".format(delta_l))
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
        startf = self.f1-4000.0
        stopf  = self.f1+4000.0
        if _kat.maxtem == "off":
            nmStr = None
        else:
            nmStr = "0 0"
        code = """
ad f1p {0} {1} nPRM2
ad f1m {0} -{1} nPRM2
xaxis mod1 f lin {2} {3} 200
put f1p f $x1 
put f1m f $mx1 
""".format(nmStr, self.f1, startf, stopf)
        kat.parseCommands(code)
        out = kat.run()
        ax.plot(out.x-self.f1,np.abs(out["f1p"]), label=" f1")
        ax.plot(out.x-self.f1,np.abs(out["f1m"]), label="-f1")
        ax.set_xlim([np.min(out.x-self.f1), np.max(out.x-self.f1)])
        ax.set_xlabel("delta_f1 [Hz]")
        ax.set_ylabel('sqrt(W) ')
        ax.grid(True)
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
        self.compute_derived_resonances(kat)

    def compute_derived_resonances(self, kat):
        self.fsrX = 0.5 * clight / float(kat.LX.L)
        self.fsrY = 0.5 * clight / float(kat.LY.L)
        self.fsrPRC = 0.5 * clight / self.lPRC
        self.fsrSRC = 0.5 * clight / self.lSRC
        self.f1_PRC = 3.5 * self.fsrPRC
        
    def lengths_status(self, kat):
        print(" .--------------------------------------------------.")
        print("| - arm length:                                     |")
        print("| Lx   = {:11.7}m, Ly   = {:11.7}m          |".format(float(kat.LX.L), float(kat.LY.L)))
        print("| - small MI and recycling lengths:                 | ")
        print("| lx   = {:11.7}m, ly   = {:11.7}m          |".format(self.lx, self.ly))
        print("| lpr  = {:11.7}m, lsr  = {:11.7}m          |".format(self.lpr, self.lsr))
        print("| lMI  = {:11.7}m, lSchnupp = {:11.5}m      |".format(self.lMI, self.lSchnupp))
        print("| lPRC = {:11.7}m, lSRC = {:11.7}m          |".format(self.lPRC, self.lSRC))
        print("+---------------------------------------------------+")
        print("| - associated cavity frequencies [Hz]:             |")
        print("| fsrx   = {:11.8}, fsry   = {:11.8}        |".format(self.fsrX, self.fsrY))
        print("| fsrPRC = {:11.8}, fsrSRC = {:11.8}        |".format(self.fsrPRC, self.fsrSRC))
        print("| f1_PRC = {:11.8}                              |".format(self.f1_PRC))
        print("| f1     = {:11.8}, f2     = {:11.9}        |".format(self.f1, self.f2))
        print(" `--------------------------------------------------'")
                    
    def apply_lock_feedback(self, kat, out):
        tuning = get_tunings(kat, self.tunings)
        if "ETMX_lock" in out.ylabels:
            tuning["ETMX"] += float(out["ETMX_lock"])
        else:
            print(" ** Warning: could not find ETMX lock")
        if "ETMY_lock" in out.ylabels:
            tuning["ETMY"] += float(out["ETMY_lock"])
        else:
            print(" ** Warning: could not find ETMY lock")
        if "PRCL_lock" in out.ylabels:
            tuning["PRM"]  += float(out["PRCL_lock"])
        else:
            print(" ** Warning: could not find PRCL lock")
        if ("MICH_lock" in out.ylabels) and ("ITMY_lock" in out.ylabels):
            tuning["ITMX"] += float(out["MICH_lock"])
            tuning["ITMY"] += float(out["ITMY_lock"])
        else:
            print(" ** Warning: could not find MICH (ITMY) lock")
        if "SRCL_lock" in out.ylabels:
            tuning["SRM"]  += float(out["SRCL_lock"])
        else:
            print(" ** Warning: could not find SRCL lock")
        set_tunings(kat, tuning)

        
    def pretune(self, _kat, pretune_precision=1.0e-4, verbose=False):
        print("-- pretuning interferometer to precision {0:2g} deg = {1:2g} m".format(pretune_precision, pretune_precision*_kat.lambda0/360.0))
        kat=_kat.deepcopy()
        #_kat.ITMX.phi=0.0
        #_kat.ITMY.phi=0.0
        vprint(verbose, "   scanning X arm (maximising power)")
        kat1 = kat.deepcopy()
        make_transparent(kat1,["PRM","SRM"])
        make_transparent(kat1,["ITMY", "ETMY"])
        kat1.BS.setRTL(0.0,1.0,0.0) # set BS refl. for X arm
        phi, precision = self.scan_to_precision(kat1, self.preARMX, pretune_precision)
        phi=round(phi/pretune_precision)*pretune_precision
        phi=round_to_n(phi,5)
        vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
        self.preARMX.apply_tuning(_kat,phi)
    
        vprint(verbose, "   scanning Y arm (maximising power)")
        kat = _kat.deepcopy()
        make_transparent(kat,["PRM","SRM"])
        make_transparent(kat,["ITMX", "ETMX"])
        kat.BS.setRTL(1.0,0.0,0.0) # set BS refl. for Y arm
        phi, precision = self.scan_to_precision(kat, self.preARMY, pretune_precision)
        phi=round(phi/pretune_precision)*pretune_precision
        phi=round_to_n(phi,5)
        vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
        self.preARMY.apply_tuning(_kat,phi)
    
        vprint(verbose, "   scanning MICH (minimising power)")
        kat = _kat.deepcopy()
        make_transparent(kat,["PRM","SRM"])
        phi, precision = self.scan_to_precision(kat, self.preMICH, pretune_precision, minmax="min", precision=30.0)
        phi=round(phi/pretune_precision)*pretune_precision
        phi=round_to_n(phi,5)
        vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
        self.preMICH.apply_tuning(_kat,phi, add=True)

        vprint(verbose, "   scanning PRCL (maximising power)")
        kat = _kat.deepcopy()
        make_transparent(kat,["SRM"])
        phi, precision = self.scan_to_precision(kat, self.prePRCL, pretune_precision)
        phi=round(phi/pretune_precision)*pretune_precision
        phi=round_to_n(phi,5)
        vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
        self.prePRCL.apply_tuning(_kat,phi)

        vprint(verbose, "   scanning SRCL (maximising carrier power, then adding 90 deg)")
        kat = _kat.deepcopy()
        phi, precision = self.scan_to_precision(kat, self.preSRCL, pretune_precision, phi=0)
        phi=round(phi/pretune_precision)*pretune_precision
        phi=round_to_n(phi,4)-90.0
        vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
        self.preSRCL.apply_tuning(_kat,phi)
        print("   ... done")
        
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
    
        tunings = get_tunings(kat, self.tunings)
        _maxtemStr = "{:3}".format(tunings["maxtem"])
        if tunings["maxtem"] == -1:
            _maxtemStr="off"
        print(" .--------------------------------------------------.")
        print(" | pretuned for maxtem = {}, phase = {:2}            |".format(_maxtemStr, int(kat.phase)))
        keys_t = list(tunings.keys())
        keys_t.remove("maxtem")
        print(" .--------------------------------------------------.")
        print(" | port   power[W] pow. ratio | optics   tunings    |")
        print(" +----------------------------|---------------------+")
        idx_p = 0
        idx_t = 0
        while (idx_p < len(pretune_DOFs) or idx_t < len(keys_t)):
            if idx_p < len(pretune_DOFs):
                p = pretune_DOFs[idx_p]
                print(" | {:5}: {:9.4g} {:9.4g} |".format(p.name, float(out[p.port.name]), float(out[p.port.name])/Pin),end="")
                idx_p +=1
            else:
                print(" |                            |", end="")
            if idx_t < len(keys_t):
                t=keys_t[idx_t]
                print(" {:5}: {:9.3g}    |".format(t, float(self.tunings[t])))
                idx_t +=1
            else:
                print("                     |")
        print(" `--------------------------------------------------'")

    # probably extra and can be removed
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
            ax.grid(True)
        plt.tight_layout()
        plt.show(block=0)

    def _strToDOFs(self, DOFs):
        dofs = []
        
        for _ in DOFs:
            if isinstance(_, six.string_types):
                if _ in self.__DOFs:
                    dofs.append(self.__DOFs[_])
                else:
                    raise pkex.BasePyKatException("Could not find DOF called `%s`. Possible DOF options: %s" % (_, str(list(self.__DOFs.keys()))))
            else:
                raise pkex.BasePyKatException("'%s' not possible DOF options: %s" % (_, str(list(self.__DOFs.keys()))))
        
        return dofs

    def plot_error_signals(self, _kat, xlimits=[-1,1], DOFs=None, plotDOFs=None,
                                replaceDOFSignals=False, block=0, fig=None, legend=None):
        """
        Displays error signals for a given kat file. Can also be used to plot multiple
        DOF's error signals against each other for visualising any cross coupling.
        
        _kat: LIGO-like kat object.
        xlimits: Range of DOF to plot in degrees
        DOFs: list, DOF names to compute. Default: DARM, CARM, PRCL, SRCL, MICH
        plotDOFs: list, DOF names to plot against each DOF. If None the same DOF as in DOFs is plotted.
        block: Boolean, for plot blocking terminal or not if being shown
        replaceDOFSignals: Bool, replaces already present signals for any DOF if already defined in kat. Regardless of this value, it will add default signals if none found.
        fig: figure, uses predefined figure, when defined it won't be shown automatically
        legend: string, if no plotDOFs is defined this legend is shown
        
        Example:
            import pykat
            from pykat.gw_detectors import ifo

            ligo = ifo.aLIGO()
            
            # Plot default
            ligo.plot_error_signals(ligo.kat, block=True)
            # plot MICH and CARM against themselves
            ligo.plot_error_signals(ligo.kat, DOFs=["MICH", "CARM"], block=True)
            # plot DARM and CARM against MICH
            ligo.plot_error_signals(ligo.kat, DOFs=["MICH"], plotDOFs=["DARM", "CARM"], block=True)
        """
        
        kat = _kat.deepcopy()
        kat.verbose = False
        kat.noxaxis = True
        
        if DOFs is None:
            dofs = [self.DARM, self.CARM, self.PRCL, self.SRCL, self.MICH]
        else:
            dofs = self._strToDOFs(DOFs)
        
        # add in signals for those DOF to plot
        for _ in dofs:
            if not (not replaceDOFSignals and hasattr(kat, _.signal_name(kat))):
                kat.parseCommands(_.signal(kat))
                
        toShow = None
        
        if plotDOFs is not None:
            toShow = self._strToDOFs(plotDOFs)
        
            # Check if other DOF signals we need to include for plotting
            for _ in toShow:
                if not (not replaceDOFSignals and hasattr(kat, _.signal_name(kat))):
                    kat.parseCommands(_.signal(kat))
                    
        if fig is not None:
            _fig = fig
        else:
            _fig = plt.figure()
        
        nrows = 2
        ncols = 3
        
        if DOFs is not None:
            n = len(DOFs)
            
            if n < 3:
                nrows = 1
                ncols = n
        
        for d, idx in zip(dofs, range(1, len(dofs)+1)):
            ax = _fig.add_subplot(nrows, ncols, idx)
            
            scan_cmd = scan_optics_string(d.optics, d.factors, "scan", linlog="lin",
                                            xlimits=np.multiply(d.scale, xlimits), steps=200,
                                            axis=1, relative=True)
            kat.parseCommands(scan_cmd)
            out = kat.run()
            
            if toShow is None:
                ax.plot(out.x, out[d.signal_name(kat)], label=legend)
            else:
                for _ in toShow:
                    if legend is None:
                        legend = _.name
                        
                    ax.plot(out.x, out[_.signal_name(kat)], label=legend)
                
            ax.set_xlim([np.min(out.x), np.max(out.x)])
            ax.set_xlabel("{} [deg]".format(d.name))
            
            if plotDOFs is None:
                ax.set_ylabel('{} [W] '.format(d.port.name))
            else:
                ax.set_ylabel('Error signal [W]')
            
            ax.grid(True)
        
        if toShow is not None or legend is not None:
            plt.legend(loc=0)
           
        plt.tight_layout()
        
        if fig is None:
            plt.show(block=block)

    def set_DC_offset(self, _kat, DCoffset=None, verbose=False):
        if DCoffset:
            self.DCoffset=DCoffset
            print("-- applying user-defined DC offset:")
            pretuning = get_tunings(_kat, self.tunings)
            pretuning["ETMY"] += self.DCoffset
            pretuning["ETMX"] -= self.DCoffset
            set_tunings(_kat, pretuning)        
            kat = _kat.deepcopy()
            sigStr = self.AS_DC.signal(kat)
            signame = self.AS_DC.signal_name(kat)
            kat.parseCommands(sigStr)
            kat.noxaxis=True
            out = kat.run()
            self.DCoffsetW=float(out[signame])
        else:
            # Finding light power in AS port (mostly due to RF sidebands now
            kat = _kat.deepcopy()
            sigStr = self.AS_DC.signal(kat)
            signame = self.AS_DC.signal_name(kat)
            kat.parseCommands(sigStr)
            kat.noxaxis=True
            out = kat.run()
            print("-- adjusting DCoffset based on light in dark port:")
            waste_light = round(float(out[signame]),1)
            print("   waste light in AS port of {:2} W".format(waste_light))
            kat_lock = _kat.deepcopy()
            self.find_DC_offset(kat_lock, 2*waste_light)
            pretuning = get_tunings(kat_lock, self.tunings)
            pretuning["ETMY"] += self.DC_offset
            pretuning["ETMX"] -= self.DC_offset
            set_tunings(_kat, pretuning)
        self.DCoffset_meter = self.DCoffset / 360.0 * _kat.lambda0 
        vprint(verbose, "   DCoffset = {:6.4} deg ({:6.4}m)".format(self.DCoffset, self.DCoffset_meter))
        vprint(verbose, "   at dark port power: {:6.4}W".format(self.DCoffsetW))
    

    def find_DC_offset(self, _kat, AS_power, precision=1e-4, verbose=False):
        """
        Returns the DC offset of DARM that corrponds to the
        specified power in the AS power.
        """
        vprint(verbose, "   finding DC offset for AS power of {:3g} W".format(AS_power))
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

        vprint(verbose, "   starting peak search...")
        out=fmin(powerDiff,0,xtol=precision,ftol=1e-3,args=(kat, Xphi, Yphi, AS_power), disp=verbose)
        vprint(verbose, "   ... done")
        vprint(verbose, "   DC offset for AS_DC={} W is: {}".format(AS_power, out[0]))
        self.DCoffset = round(out[0],6)
        self.DCoffsetW = AS_power
        return self.DCoffset

    def generate_tuning_block(self, kat):
        code1 = """###########################################################################
const phi_ITMX {:.8}
const phi_ITMY {:.8}
const phi_ETMX {:.8}
const phi_ETMY {:.8}""".format(float(kat.ITMX.phi), float(kat.ITMY.phi), float(kat.ETMX.phi),float(kat.ETMY.phi))

        code2 = """
const phi_BS   {:.8}
const phi_PRM  {:.8}
const phi_SRM  {:.8}
maxtem {}
phase {}
###########################################################################""".format(float(kat.BS.phi), float(kat.PRM.phi), float(kat.SRM.phi), kat.maxtem, kat.phase)

        return "".join([code1, code2])

    def generate_errsig_block(self, kat, noplot=False):
        code1 = "###########################################################################\n" 
        sigDARM = self.DARM.signal(kat)
        sigCARM = self.CARM.signal(kat)
        sigPRCL = self.PRCL.signal(kat)
        sigMICH = self.MICH.signal(kat)
        sigSRCL = self.SRCL.signal(kat)
        code2 = "\n".join([sigDARM, sigCARM, sigPRCL, sigMICH, sigSRCL])

        code3= ""
        if noplot:
            nameDARM = self.DARM.signal_name(kat)
            nameCARM = self.CARM.signal_name(kat)
            namePRCL = self.PRCL.signal_name(kat)
            nameMICH = self.MICH.signal_name(kat)
            nameSRCL = self.SRCL.signal_name(kat)
            code3 = """
noplot {}
noplot {}
noplot {}
noplot {}
noplot {}""".format(nameDARM, nameCARM, namePRCL, nameMICH, nameSRCL)
        code4="\n###########################################################################"
        return "".join([code1, code2, code3, code4])

            
        
    def generate_locks(self, _kat, gainsAdjustment = [0.5, 0.005, 1.0, 0.5, 0.025], gains=None, accuracies=None, verbose=True):
        """
        gainsAdjustment: factors to apply to loop gains computed from optical gains
        gains: override loop gain [W per deg]
        accuracies: overwrite error signal threshold [W]
        """
        kat = _kat.deepcopy()
        # optical gains in W/rad
        ogDARM = optical_gain(kat, self.DARM, self.DARM)
        ogCARM = optical_gain(kat, self.CARM, self.CARM)
        ogPRCL = optical_gain(kat, self.PRCL, self.PRCL)
        ogMICH = optical_gain(kat, self.MICH, self.MICH)
        ogSRCL = optical_gain(kat, self.SRCL, self.SRCL)

        if gains == None:            
            # manually tuning relative gains
            factor = -1.0 * 180 / math.pi # convert from rad/W to -1 * deg/W
            gainDARM = round_to_n(gainsAdjustment[0] * factor / ogDARM, 2) # manually tuned
            gainCARM = round_to_n(gainsAdjustment[1] * factor / ogCARM, 2) # factor 0.005 for better gain hirarchy with DARM
            gainPRCL = round_to_n(gainsAdjustment[2] * factor / ogPRCL, 2) # manually tuned
            gainMICH = round_to_n(gainsAdjustment[3] * factor / ogMICH, 2) # manually tuned
            gainSRCL = round_to_n(gainsAdjustment[4] * factor / ogSRCL, 2) # gain hirarchy with MICH
            gains = [ gainDARM, gainCARM, gainPRCL, gainMICH, gainSRCL]
        self.lockGains = copy.deepcopy(gains)
        
        # rms: loop accuracies in meters (manually tuned for the loops to work
        # with the default file)
        # to compute accuracies from rms, we convert
        # rms to radians as rms_rad = rms * 2 pi/lambda
        # and then multiply by the optical gain.
        rms = [1e-13, 1e-13, 1e-12, 1e-11, 50e-11] # default accuracies in meters
        factor = 2.0 * math.pi / kat.lambda0 # convert from m to radians
        if accuracies == None:
            accDARM = round_to_n(np.abs(factor * rms[0] * ogDARM),2) 
            accCARM = round_to_n(np.abs(factor * rms[1] * ogCARM),2) 
            accPRCL = round_to_n(np.abs(factor * rms[2] * ogPRCL),2) 
            accMICH = round_to_n(np.abs(factor * rms[3] * ogMICH),2)
            accSRCL = round_to_n(np.abs(factor * rms[4] * ogSRCL),2) 
            accuracies = [accDARM, accCARM, accPRCL, accMICH, accSRCL]
        self.lockAccuracies = copy.deepcopy(accuracies)


        nameDARM = self.DARM.signal_name(kat)
        nameCARM = self.CARM.signal_name(kat)
        namePRCL = self.PRCL.signal_name(kat)
        nameMICH = self.MICH.signal_name(kat)
        nameSRCL = self.SRCL.signal_name(kat)
        self.lockNames = [nameDARM, nameCARM, namePRCL, nameMICH, nameSRCL]
        factor1 = 2.0 * math.pi / 360.0 
        factor2 = 2.0 * math.pi / kat.lambda0 
        factor3 = 360.0  / kat.lambda0
        factor4 = -1.0 * 180 / math.pi 

        if verbose:
            print(" .--------------------------------------------------.")
            print(" | Parameters for locks:                            |")
            print(" +--------------------------------------------------+")
            print(" | -- optical gains [W/rad], [W/deg] and [W/m]:     |")
            print(" | DARM: {:12.5}, {:12.5}, {:12.5}   |".format(ogDARM, ogDARM*factor1, ogDARM*factor2))
            print(" | CARM: {:12.5}, {:12.5}, {:12.5}   |".format(ogCARM, ogCARM*factor1, ogCARM*factor2))
            print(" | PRCL: {:12.5}, {:12.5}, {:12.5}   |".format(ogPRCL, ogPRCL*factor1, ogPRCL*factor2))
            print(" | MICH: {:12.5}, {:12.5}, {:12.5}   |".format(ogMICH, ogMICH*factor1, ogMICH*factor2))
            print(" | SRCL: {:12.5}, {:12.5}, {:12.5}   |".format(ogSRCL, ogSRCL*factor1, ogSRCL*factor2))
            print(" +--------------------------------------------------+")
            print(" | -- defult loop accuracies [deg], [m] and [W]:    |")
            print(" | DARM: {:12.6}, {:12.6}, {:12.6}   |".format(factor3*rms[0], rms[0], np.abs(rms[0]*ogDARM*factor2)))
            print(" | CARM: {:12.6}, {:12.6}, {:12.6}   |".format(factor3*rms[1], rms[1], np.abs(rms[1]*ogCARM*factor2)))
            print(" | PRCL: {:12.6}, {:12.6}, {:12.6}   |".format(factor3*rms[2], rms[2], np.abs(rms[2]*ogPRCL*factor2)))
            print(" | MICH: {:12.6}, {:12.6}, {:12.6}   |".format(factor3*rms[3], rms[3], np.abs(rms[3]*ogMICH*factor2)))
            print(" | SRCL: {:12.6}, {:12.6}, {:12.6}   |".format(factor3*rms[4], rms[4], np.abs(rms[4]*ogSRCL*factor2)))
            print(" +--------------------------------------------------+")
            print(" | -- extra gain factors (factor * 1/optical_gain): |")
            print(" | DARM: {:5.4} * {:12.6} = {:12.6}        |".format(gainsAdjustment[0],factor4/ogDARM, gainsAdjustment[0]*factor4/ogDARM))
            print(" | CARM: {:5.4} * {:12.6} = {:12.6}        |".format(gainsAdjustment[1],factor4/ogCARM, gainsAdjustment[1]*factor4/ogCARM))
            print(" | PRCL: {:5.4} * {:12.6} = {:12.6}        |".format(gainsAdjustment[2],factor4/ogPRCL, gainsAdjustment[2]*factor4/ogPRCL))
            print(" | MICH: {:5.4} * {:12.6} = {:12.6}        |".format(gainsAdjustment[3],factor4/ogMICH, gainsAdjustment[3]*factor4/ogMICH))
            print(" | SRCL: {:5.4} * {:12.6} = {:12.6}        |".format(gainsAdjustment[4],factor4/ogSRCL, gainsAdjustment[4]*factor4/ogSRCL))
            print(" `--------------------------------------------------'")


        
    def generate_lock_block(self, kat, verbose=False):
        if self.lockNames == None or self.lockAccuracies == None or self.lockGains == None:
            raise pkex.BasePyKatException("run generate_locks before generate_lock_block")            
        code1 = """###########################################################################
set AS_f2_I_re {} re
set CARM_err {} re
set PRCL_err {} re
set MICH_err {} re
set SRCL_err {} re
func DARM_err = $AS_f2_I_re - {}
""".format(self.lockNames[0], self.lockNames[1], self.lockNames[2], self.lockNames[3], self.lockNames[4], self.DCoffsetW)
       
        code2 = """
lock DARM_lock $DARM_err {:8.2} {:8.2}
lock CARM_lock $CARM_err {:8.2g} {:8.2g} 
lock PRCL_lock $PRCL_err {:8.2g} {:8.2g}
lock MICH_lock $MICH_err {:8.2g} {:8.2g} 
lock SRCL_lock $SRCL_err {:8.2g} {:8.2g}""".format(self.lockGains[0], self.lockAccuracies[0], self.lockGains[1], self.lockAccuracies[1], self.lockGains[2], self.lockAccuracies[2], self.lockGains[3], self.lockAccuracies[3], self.lockGains[4], self.lockAccuracies[4])
        
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
"""
        if verbose:
            print(" .--------------------------------------------------.")
            print(" | Lock commands used:                              |")
            print(" +--------------------------------------------------+")
            for l in code2.splitlines():
                print (" | {:49}|".format(l))
            print(" `--------------------------------------------------'")
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
                                
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

def get_tunings(kat, tunings):
    """
    returns the tunings of optical components and the corresponding values
    for maxtem and phase
    """
    keys = list(tunings.keys())
    if "maxtem" in keys:
        tunings["maxtem"] = kat.maxtem
        keys.remove("maxtem")
    if "phase" in keys:
        tunings["phase"] = kat.phase
        keys.remove("phase")
    for comp in keys:
        tunings[comp] = kat.components[comp].phi
    return tunings
    
def set_tunings(kat, tunings):
    """
    sets the tunings of optical components and the corresponding values
    for maxtem and phase
    """
    keys = list(tunings.keys())
    if "maxtem" in keys:
        kat.maxtem = tunings["maxtem"]
        keys.remove("maxtem")
    if "phase" in keys:
        kat.phase = tunings["phase"] 
        keys.remove("phase")
    for comp in keys:
        kat.components[comp].phi = tunings[comp] 

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

def reconnect_nodes(kat, component1, idx1, node_name, verbose=False):
    c_string = component1.getFinesseText()
    c = c_string[0].split()
    new_string = " ".join(c[:-2])
    nodes = {}
    nodes[0] = c[-2]
    nodes[1] = c[-1]
    nodes[idx1]=node_name
    new_string = new_string + " " + nodes[0] + " " + nodes[1]
    vprint(verbose, "   new string ='{}'".format(new_string))
    kat.parseCommands(new_string)

def remove_commands(kat, _commands, verbose=False):
    commands=make_list_copy(_commands)
    # removing commands
    for o in kat.commands.values():
        if o.name in commands:
            o.remove()
            commands = [c for c in commands if c != o.name]
            vprint(verbose, '   {} removed'.format(o))
    if len(commands) != 0:
        raise pkex.BasePyKatException("Cannot find command(s) {}".format(commands))
    return kat
    
def remove_components(kat, _components, component_in=None, component_out=None, verbose=False):
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
            vprint(verbose, '  {} removed'.format(o))
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
    
    scan_string = scan_optics_string(DOF.optics, DOF.factors, "scan", linlog="lin",
                                     xlimits=xlimits, steps=steps, axis=1, relative=relative)
    
    kat.parseCommands(scan_string)
    kat.parseCommands(DOF.signal(kat))
    
    return kat.run()

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
    
    scan_string = scan_optics_string(optics, factors, "scan", linlog="lin", xlimits=xlimits,
                                      steps=steps, axis=1,relative=relative)
    
    kat.parseCommands(scan_string)
    
    return kat.run()

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
    return float(np.real(out[_detName]))

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
    if debug:
        print("  stepsize (precision) of scan: {0:g}".format(stepsize))

    _max, _min = peak.peakdetect( out[detector],out.x, 1)
    
    if debug:
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
        
    if debug:
        plt.plot(X_out,Y_out,'o')
        plt.xlabel('tuning [deg]')
        plt.ylabel('{0} output'.format(detector))
        plt.show(block=0)
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


def vprint(verbose, printstr):
    if verbose:
        print(printstr)
