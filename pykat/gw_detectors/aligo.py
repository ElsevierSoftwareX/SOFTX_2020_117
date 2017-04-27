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
from . import IFO, DOF, Port, vprint, clight, nsilica, find_peak, make_transparent, reconnect_nodes, remove_commands, remove_components, round_to_n, scan_optics_string

import pykat.components
import pykat.exceptions as pkex
import pykat.external.peakdetect as peak
import matplotlib.pyplot as plt
import pkg_resources

from scipy.optimize import fmin

class ALIGO_IFO(IFO):   
    """
    This contains aLIGO specific methods for computing interferometer
    variables.
    
    Functions that operate on the kat/IFO objects, manipulate their
    structure and return a new one, should not be included here. They should
    be separate functions that are called by the user. Functions that are
    chained together by the user, like pretuning, should be separate functions.
    
    The functions here should be those that update the kat with information
    from the IFO or vice-versa.
    """
    def compute_derived_resonances(self):
        self.fsrX = 0.5 * clight / float(self.kat.LX.L)
        self.fsrY = 0.5 * clight / float(self.kat.LY.L)
        self.fsrPRC = 0.5 * clight / self.lPRC
        self.fsrSRC = 0.5 * clight / self.lSRC
        self.f1_PRC = 3.5 * self.fsrPRC
    
    def compute_derived_lengths(self, verbose=False):
        """
        Compute derived length from individual space components.
        Design values are currently:
        lPRC = 57.656, lSRC = 56.008, lSchnupp = 0.08
        and the individual lengths:
        PRC: Lp1 16.6107, Lp2 16.1647, Lp3 19.5381 
        SRC: Ls1 15.7586, Ls2 15.4435, Ls3 19.3661
        """
        # distances between HR surfaces:
        self.lpr = self.kat.lp1.L + self.kat.lp2.L + self.kat.lp3.L
        self.lx = self.kat.lx1.L + self.kat.BSsub1.L * self.kat.BSsub1.n + self.kat.ITMXsub.L * self.kat.ITMXsub.n
        self.ly = self.kat.ly1.L + self.kat.ITMYsub.L * self.kat.ITMYsub.n
        self.lsr = self.kat.ls1.L + self.kat.ls2.L + self.kat.ls3.L + self.kat.BSsub2.L * self.kat.BSsub2.n
        # resulting combined distances (single, not roundtrip)
        self.lMI =  0.5 * (self.lx + self.ly)
        self.lPRC = self.lpr + self.lMI
        self.lSRC = self.lsr + self.lMI
        self.lSchnupp = self.lx - self.ly
    
        self.compute_derived_resonances()
        
    def lengths_status(self):
        self.compute_derived_lengths()
        
        print(" .--------------------------------------------------.")
        print("| - arm length:                                     |")
        print("| Lx   = {:11.7}m, Ly   = {:11.7}m          |".format(float(self.kat.LX.L), float(self.kat.LY.L)))
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
    
    def adjust_PRC_length(self, verbose=False):
        """
        Adjust PRC length so that it fulfils the requirement
        lPRC = (N+1/2) * c/(2*f1), see [1] equation C.1
        In the current design N=3.
    
        Function directly alters the associated kat object.
        """
        kat = self.kat
        
        vprint(kat.verbose, "-- adjusting PRC length")
        ltmp = 0.5 * clight / kat.IFO.f1
        delta_l = 3.5 * ltmp - kat.IFO.lPRC
        vprint(kat.verbose, "   adusting kat.lp1.L by {:.4g}m".format(delta_l))
        kat.lp1.L += delta_l
    
        kat.IFO.compute_derived_lengths(kat)
        
 
def assert_aligo_ifo_kat(kat):
    if not isinstance(kat.IFO, ALIGO_IFO):
        raise pkex.BasePyKatException("\033[91mkat file is not an ALIGO_IFO compatiable kat\033[0m")
              
def make_kat(name="default", katfile=None, verbose = False, debug=False):
    """
    Returns a kat object and fills in the kat.IFO property for storing
    the associated interferometer data.
    """
    
    names = ['default', 'LLO', 'LHO']
    
    if debug:
        kat = finesse.kat(tempdir=".",tempname="test")
    else:
        kat = finesse.kat()
        
    kat.verbose=verbose
    
    # Create empty object to just store whatever DOFs, port, variables in
    # that will be used by processing functions
    kat.IFO = ALIGO_IFO(kat,
                        # Define which keys are used for a tuning description
                        ["maxtem", "phase"],
                        # Define which mirrors create the tuning description
                        ["PRM", "ITMX", "ETMX", "ITMY", "ETMY", "BS", "SRM"])
    
    kat.IFO._data_path=pkg_resources.resource_filename('pykat.gw_detectors', 'finesse_files/')

    kat.IFO.rawBlocks = BlockedKatFile()
    
    if katfile:
        kat.load(katfile)
        kat.IFO.rawBlocks.read(katfile)
    else:
        """
        if name not in names: # TODO different files not yet implemented
            printf("aLIGO name `{}' not recognised, must be 'default', 'LLO' or 'LHO'",name)
        """
        if name != "default":
            printf("aLIGO name `{}' not recognised, using 'default'", name)
        
        kat.load(kat.IFO._data_path+"aLIGO.kat")
        kat.IFO.rawBlocks.read(kat.IFO._data_path+"aLIGO.kat")
        
    # ----------------------------------------------------------------------
    # set variables to zero first
    kat.IFO.DCoffset = 0.0
    kat.IFO.DCoffsetW = 0.0
    
    # ----------------------------------------------------------------------
    # get and derive parameters from the kat file
    
    # get main frequencies
    if "f1" in kat.constants.keys():
        kat.IFO.f1 = float(kat.constants["f1"].value)
    else:
        kat.IFO.f1 = 9099471.0
        
    if "f2" in kat.constants.keys():
        kat.IFO.f2 = float(kat.constants["f2"].value)
    else:
        kat.IFO.f2 = 5.0 * kat.IFO.f1
        
    if "f3" in kat.constants.keys():
        kat.IFO.f3 = float(kat.constants["f3"].value)
        
    # TODO add else here!
    # check modultion frequencies
    if (5 * kat.IFO.f1 != kat.IFO.f2):
        print(" ** Warning: modulation frequencies do not match: 5*f1!=f2")
    
    # defining a dicotionary for the main mirror positions (tunings),
    # keys should include maxtem, phase and all main optics names
    #kat.IFO.tunings = get_tunings(dict.fromkeys(["maxtem", "phase", "PRM", "ITMX", "ETMX", "ITMY", "ETMY", "BS", "SRM"]))
    kat.IFO.compute_derived_lengths()
        
    # ----------------------------------------------------------------------
    # define ports and signals 
    
    # useful ports
    kat.IFO.POP_f1  = Port(kat.IFO, "POP_f1",  "nPOP",  kat.IFO.f1, phase=101)
    kat.IFO.POP_f2  = Port(kat.IFO, "POP_f2",  "nPOP",  kat.IFO.f2, phase=13)
    kat.IFO.REFL_f1 = Port(kat.IFO, "REFL_f1", "nREFL", kat.IFO.f1, phase=101)
    kat.IFO.REFL_f2 = Port(kat.IFO, "REFL_f2", "nREFL", kat.IFO.f2, phase=14)
    kat.IFO.AS_DC   = Port(kat.IFO, "AS_DC", "nSRM2")
    kat.IFO.POW_BS  = Port(kat.IFO, "PowBS", "nPRBS*")
    kat.IFO.POW_X   = Port(kat.IFO, "PowX",  "nITMX2")
    kat.IFO.POW_Y   = Port(kat.IFO, "PowY",  "nITMY2")

    # pretune DOF
    kat.IFO.preARMX =  DOF(kat.IFO, "ARMX", kat.IFO.POW_X,   "", "ETMX", 1, 1.0)
    kat.IFO.preARMY =  DOF(kat.IFO, "ARMY", kat.IFO.POW_Y,   "", "ETMY", 1, 1.0)
    kat.IFO.preMICH =  DOF(kat.IFO, "AS"  , kat.IFO.AS_DC,   "", ["ITMX", "ETMX", "ITMY", "ETMY"], [1,1,-1,-1], 6.0)
    kat.IFO.prePRCL =  DOF(kat.IFO, "PRCL", kat.IFO.POW_BS,  "", "PRM",  1, 10.0)
    kat.IFO.preSRCL =  DOF(kat.IFO, "SRCL", kat.IFO.AS_DC,   "", "SRM",  1, 10.0)
    
    # control scheme as in [1] Table C.1  
    kat.IFO.PRCL =  DOF(kat.IFO, "PRCL", kat.IFO.POP_f1,  "I", "PRM", 1, 100.0)
    kat.IFO.MICH =  DOF(kat.IFO, "MICH", kat.IFO.POP_f2,  "Q", ["ITMX", "ETMX", "ITMY", "ETMY"], [1,1,-1,-1], 100.0)
    kat.IFO.CARM =  DOF(kat.IFO, "CARM", kat.IFO.REFL_f1, "I", ["ETMX", "ETMY"], [1, 1], 1.5)
    kat.IFO.DARM =  DOF(kat.IFO, "DARM", kat.IFO.AS_DC,   "",  ["ETMX", "ETMY"], [1,-1], 1.0)
    kat.IFO.SRCL =  DOF(kat.IFO, "SRCL", kat.IFO.REFL_f2, "I", "SRM", 1, 1e2)
    
    kat.IFO.DOFs = {}
    
    for _ in inspect.getmembers(kat.IFO, lambda x: isinstance(x, DOF)):
        kat.IFO.DOFs[_[0]] = _[1]
        
    kat.IFO.Ports = {}
    
    for _ in inspect.getmembers(kat.IFO, lambda x: isinstance(x, Port)):
        kat.IFO.Ports[_[0]] = _[1]

    kat.IFO.lockNames = None
    
    return kat
    


def plot_f1_PRC_resonance(_kat, ax=None, show=True):
    """
    Plot the sideband amplitudes for modulation frequecy
    f1 (~ 9MHz) in the PRC, to check the resonance
    condition.
    """
    assert_aligo_ifo_kat(_kat)
    
    kat = _kat.deepcopy()
    
    # Don't need locks for this plot so remove if present
    kat.removeBlock('locks', False)
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    startf = kat.IFO.f1 - 4000.0
    stopf  = kat.IFO.f1 + 4000.0

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
            """.format(nmStr, kat.IFO.f1, startf, stopf)
            
    kat.parseCommands(code)
    
    out = kat.run()
    
    ax.plot(out.x-kat.IFO.f1,np.abs(out["f1p"]), label=" f1")
    ax.plot(out.x-kat.IFO.f1,np.abs(out["f1m"]), label="-f1", ls="--")
    ax.set_xlim([np.min(out.x-kat.IFO.f1), np.max(out.x-kat.IFO.f1)])
    ax.set_xlabel("delta_f1 [Hz]")
    ax.set_ylabel('sqrt(W) ')
    ax.grid(True)
    ax.legend()
    ax.figure.set_tight_layout(True)
    
    if show: plt.show()
    
def apply_lock_feedback(kat, out):
    tuning = kat.IFO.get_tunings()
    
    if "ETMX_lock" in out.ylabels:
        tuning["ETMX"] += float(out["ETMX_lock"])
    else:
        pkex.printWarning(" ** Warning: could not find ETMX lock")
        
    if "ETMY_lock" in out.ylabels:
        tuning["ETMY"] += float(out["ETMY_lock"])
    else:
        pkex.printWarning(" ** Warning: could not find ETMY lock")
        
    if "PRCL_lock" in out.ylabels:
        tuning["PRM"]  += float(out["PRCL_lock"])
    else:
        pkex.printWarning(" ** Warning: could not find PRCL lock")
        
    if ("MICH_lock" in out.ylabels) and ("ITMY_lock" in out.ylabels):
        tuning["ITMX"] += float(out["MICH_lock"])
        tuning["ITMY"] += float(out["ITMY_lock"])
    else:
        pkex.printWarning(" ** Warning: could not find MICH (ITMY) lock")
        
    if "SRCL_lock" in out.ylabels:
        tuning["SRM"]  += float(out["SRCL_lock"])
    else:
        pkex.printWarning(" ** Warning: could not find SRCL lock")
        
    kat.IFO.apply_tunings(tuning)
    
def scan_to_precision(kat, DOF, pretune_precision, minmax="max", phi=0.0, precision=60.0):
    while precision > pretune_precision * DOF.scale:
        out = kat.IFO.scan_DOF(DOF, xlimits = [phi-1.5*precision, phi+1.5*precision])
        phi, precision = find_peak(out, DOF.port.portName, minmax=minmax)
        
    return phi, precision
    
    
def pretune(_kat, pretune_precision=1.0e-4, verbose=False):
    assert_aligo_ifo_kat(_kat)
    
    # This function needs to apply a bunch of pretunings to the original
    # kat and associated IFO object passed in
    IFO = _kat.IFO
    
    print("-- pretuning interferometer to precision {0:2g} deg = {1:2g} m".format(pretune_precision, pretune_precision*_kat.lambda0/360.0))
    
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    
    vprint(verbose, "   scanning X arm (maximising power)")
    
    make_transparent(kat, ["PRM", "SRM"])
    make_transparent(kat, ["ITMY", "ETMY"])
    
    kat.BS.setRTL(0.0, 1.0, 0.0) # set BS refl. for X arm
    
    phi, precision = scan_to_precision(kat, IFO.preARMX, pretune_precision)
    phi = round(phi/pretune_precision)*pretune_precision
    phi = round_to_n(phi,5)
    
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    
    IFO.preARMX.apply_tuning(phi)

    vprint(verbose, "   scanning Y arm (maximising power)")
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    
    make_transparent(kat,["PRM","SRM"])
    make_transparent(kat,["ITMX", "ETMX"])
    kat.BS.setRTL(1.0,0.0,0.0) # set BS refl. for Y arm
    phi, precision = scan_to_precision(kat, IFO.preARMY, pretune_precision)
    phi=round(phi/pretune_precision)*pretune_precision
    phi=round_to_n(phi,5)
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.preARMY.apply_tuning(phi)

    vprint(verbose, "   scanning MICH (minimising power)")
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    
    make_transparent(kat,["PRM","SRM"])
    phi, precision = scan_to_precision(kat, IFO.preMICH, pretune_precision, minmax="min", precision=30.0)
    phi=round(phi/pretune_precision)*pretune_precision
    phi=round_to_n(phi,5)
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.preMICH.apply_tuning(phi, add=True)

    vprint(verbose, "   scanning PRCL (maximising power)")
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    make_transparent(kat,["SRM"])
    phi, precision = scan_to_precision(kat, IFO.prePRCL, pretune_precision)
    phi=round(phi/pretune_precision)*pretune_precision
    phi=round_to_n(phi,5)
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.prePRCL.apply_tuning(phi)

    vprint(verbose, "   scanning SRCL (maximising carrier power, then adding 90 deg)")
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    
    phi, precision = scan_to_precision(kat, IFO.preSRCL, pretune_precision, phi=0)
    phi=round(phi/pretune_precision)*pretune_precision
    phi=round_to_n(phi,4)-90.0
    
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.preSRCL.apply_tuning(phi)
    
    print("   ... done")
    


def pretune_status(_kat):
    assert_aligo_ifo_kat(_kat)
    
    kat = _kat.deepcopy()
    kat.verbose = False
    kat.noxaxis = True
    
    pretune_DOFs = [kat.IFO.preARMX, kat.IFO.preARMY, kat.IFO.prePRCL, kat.IFO.preMICH, kat.IFO.preSRCL]
    
    _detStr=""
    
    for p in pretune_DOFs:
        _sigStr = p.port.signal(kat)
        _detStr = "\n".join([_detStr, _sigStr])
        
    kat.parseCommands(_detStr)
    out = kat.run()
    Pin = float(kat.L0.P)

    tunings = kat.IFO.get_tunings()
    
    if tunings['keys']["maxtem"] == -1:
        _maxtemStr="off"
    else:
        _maxtemStr = "{:3}".format(tunings['keys']["maxtem"])
        
    print(" .--------------------------------------------------.")
    print(" | pretuned for maxtem = {}, phase = {:2}            |".format(_maxtemStr, int(kat.phase)))
    
    keys_t = list(tunings.keys())
    keys_t.remove("keys")
    
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
            print(" {:5}: {:9.3g}    |".format(t, float(tunings[t])))
            idx_t +=1
        else:
            print("                     |")
            
    print(" `--------------------------------------------------'")

# probably extra and can be removed
def power_ratios(_kat):
    assert_aligo_ifo_kat(_kat)
    
    kat = _kat.deepcopy()
    kat.verbose = False
    kat.noxaxis = True

    ports = [kat.IFO.POW_X, kat.IFO.POW_Y, kat.IFO.AS_DC, kat.IFO.POW_BS]
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
        

def plot_pretuning_powers(self, _kat, xlimits=[-10,10]):
    assert_aligo_ifo_kat(_kat)
    
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
    
def plot_error_signals(_kat, xlimits=[-1,1], DOFs=None, plotDOFs=None,
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
        dofs = [kat.IFO.DARM, kat.IFO.CARM, kat.IFO.PRCL, kat.IFO.SRCL, kat.IFO.MICH]
    else:
        dofs = kat.IFO.strToDOFs(DOFs)
    
    # add in signals for those DOF to plot
    for _ in dofs:
        if not (not replaceDOFSignals and hasattr(kat, _.signal_name())):
            kat.parseCommands(_.signal())
            
    toShow = None
    
    if plotDOFs is not None:
        toShow = self._strToDOFs(plotDOFs)
    
        # Check if other DOF signals we need to include for plotting
        for _ in toShow:
            if not (not replaceDOFSignals and hasattr(kat, _.signal_name())):
                kat.parseCommands(_.signal())
                
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
            ax.plot(out.x, out[d.signal_name()], label=legend)
        else:
            for _ in toShow:
                if legend is None:
                    legend = _.name
                    
                ax.plot(out.x, out[_.signal_name()], label=legend)
            
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
        
        
def set_DC_offset(_kat, DCoffset=None, verbose=False):
    if DCoffset:
        _kat.IFO.DCoffset = DCoffset
        
        print("-- applying user-defined DC offset:")
        
        tunings = _kat.IFO.get_tunings()
        
        tunings["ETMY"] += _kat.IFO.DCoffset
        tunings["ETMX"] -= _kat.IFO.DCoffset
        
        _kat.IFO.apply_tunings(tunings)        
        
        kat = _kat.deepcopy()
        
        sigStr = kat.IFO.AS_DC.signal()
        signame = kat.IFO.AS_DC.signal_name()
        
        kat.parseCommands(sigStr)
        kat.noxaxis=True
        
        out = kat.run()
        
        _kat.IFO.DCoffsetW = float(out[signame])
    else:
        # Finding light power in AS port (mostly due to RF sidebands now
        kat = _kat.deepcopy()
        
        sigStr = kat.IFO.AS_DC.signal()
        signame = kat.IFO.AS_DC.signal_name()
        
        kat.parseCommands(sigStr)
        kat.noxaxis=True
        
        out = kat.run()
        
        print("-- adjusting DCoffset based on light in dark port:")
        
        waste_light = round(float(out[signame]),1)
        
        print("   waste light in AS port of {:2} W".format(waste_light))
        
        #kat_lock = _kat.deepcopy()
        
        find_DC_offset(_kat, 2*waste_light)
        
    vprint(verbose, "   DCoffset = {:6.4} deg ({:6.4}m)".format(_kat.IFO.DCoffset, _kat.IFO.DCoffset / 360.0 * _kat.lambda0 ))
    vprint(verbose, "   at dark port power: {:6.4}W".format(_kat.IFO.DCoffsetW))


def find_DC_offset(_kat, AS_power, precision=1e-4, verbose=False):
    """
    Returns the DC offset of DARM that corrponds to the
    specified power in the AS power.
    """
    vprint(verbose, "   finding DC offset for AS power of {:3g} W".format(AS_power))
    
    kat = _kat.deepcopy()
    kat.verbose = False
    kat.noxaxis = True
    
    _sigStr = kat.IFO.AS_DC.signal()
    
    kat.parseCommands(_sigStr)
    
    Xphi = float(kat.ETMX.phi)
    Yphi = float(kat.ETMY.phi)

    def powerDiff(phi, kat, Xphi, Yphi, AS_power):
        kat.ETMY.phi = Yphi + phi
        kat.ETMX.phi = Xphi - phi
        
        out = kat.run()
        
        #print(out[self.AS_DC.name]-AS_power)
        return np.abs(out[self.AS_DC.name] - AS_power)

    vprint(verbose, "   starting peak search...")
    out = fmin(powerDiff, 0, xtol=precision, ftol=1e-3, args=(kat, Xphi, Yphi, AS_power), disp=verbose)
    
    vprint(verbose, "   ... done")
    vprint(verbose, "   DC offset for AS_DC={} W is: {}".format(AS_power, out[0]))
    
    _kat.IFO.DCoffset = round(out[0],6)
    _kat.IFO.DCoffsetW = AS_power
    
    tunings = _kat.IFO.get_tunings()
    tunings["ETMY"] += _kat.IFO.DC_offset
    tunings["ETMX"] -= _kat.IFO.DC_offset
    
    _kat.IFO.apply_tunings(tunings)
    