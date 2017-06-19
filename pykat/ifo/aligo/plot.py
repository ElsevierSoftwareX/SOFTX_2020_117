import matplotlib.pyplot as plt

from . import assert_aligo_ifo_kat
from .. import scan_optics_string

import pykat.ifo
import numpy as np

def f1_PRC_resonance(_kat, ax=None, show=True):
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

def pretuning_powers(self, _kat, xlimits=[-10,10]):
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
    
def error_signals(_kat, xlimits=[-1,1], DOFs=None, plotDOFs=None,
                            replaceDOFSignals=False, block=True, fig=None, legend=None, steps=100):
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
    kat.removeBlock("locks", False)
    
    if DOFs is None:
        dofs = [kat.IFO.DARM, kat.IFO.CARM, kat.IFO.PRCL, kat.IFO.SRCL, kat.IFO.MICH]
    else:
        dofs = kat.IFO.strsToDOFs(DOFs)
    
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
        
        kat.removeBlock("SCAN", False)
        
        scan_cmd = scan_optics_string(d.optics, d.factors, "scan", linlog="lin",
                                        xlimits=np.multiply(d.scale, xlimits), steps=steps,
                                        axis=1, relative=True)
                                  
        kat.parseCommands(scan_cmd, addToBlock="SCAN")
        
        out = kat.run()
        
        if d.name == "DARM":
            DC_Offset = kat.IFO.DCoffsetW
        else:
            DC_Offset = 0
        
        if toShow is None:
            ax.plot(out.x, out[d.signal_name()] - DC_Offset, label=legend)
        else:
            for _ in toShow:
                if legend is None:
                    legend = _.name
                    
                ax.plot(out.x, out[_.signal_name()] - DC_Offset, label=legend)
            
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
        
def cavity_modes(kat, node, ax=None, show=True):
    mmx, mmy, qs = pykat.ifo.mismatch_cavities(kat, node)
    
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    
    for c, qx, qy in qs:
        ax.scatter(qx.z, qx.zr, label=c+" (x)")
        ax.scatter(qy.z, qy.zr, marker='+', label=c+" (y)")
    
    ax.set_xlabel('$z$')
    ax.set_ylabel('$z_r$')
    ax.legend(loc=0, fontsize=10, bbox_to_anchor=(1.04,1), loc="upper left")
    
    plt.tight_layout()
    
    if show: plt.show()