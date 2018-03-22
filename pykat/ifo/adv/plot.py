import matplotlib.pyplot as plt

from . import assert_adv_ifo_kat
from .. import scan_optics_string

from pykat.ifo.plot import *

import pykat.ifo
import numpy as np
import six 

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
            
    kat.parse(code)
    
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
    
def error_signals(_kat, xlimits=[-1,1], DOFs=None, plotDOFs=None, replaceDOFSignals=False,
                  block=True, fig=None, legend=None, steps=100, figsize=None):
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
        dofs = kat.IFO.LSC_DOFs
    else:
        dofs = kat.IFO.strsToDOFs(DOFs)
    
    # add in signals for those DOF to plot
    for _ in dofs:
        if not (not replaceDOFSignals and hasattr(kat, _.signal_name())):
            kat.parse(_.signal())
            
    toShow = None
    
    if plotDOFs is not None:
        toShow = self._strToDOFs(plotDOFs)
    
        # Check if other DOF signals we need to include for plotting
        for _ in toShow:
            if not (not replaceDOFSignals and hasattr(kat, _.signal_name())):
                kat.parse(_.signal())
                
    if fig is not None:
        _fig = fig
    else:
        _fig = plt.figure(figsize=figsize)
    
    n = len(dofs)
    if n < 4:
        nrows = 1
        ncols = n
    elif n == 4:
        nrows = 2
        ncols = 2
    else:
        ncols = 3
        nrows = int(np.ceil(n/ncols))
    
    for d, idx in zip(dofs, range(1, len(dofs)+1)):
        ax = _fig.add_subplot(nrows, ncols, idx)
        
        kat.removeBlock("SCAN", False)
        
        scan_cmd = scan_optics_string(d.optics, d.factors, "scan", linlog="lin",
                                        xlimits=np.multiply(d.scale, xlimits), steps=steps,
                                        axis=1, relative=True)
                                  
        kat.parse(scan_cmd, addToBlock="SCAN")
                
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
    
    
def amps_vs_dof(kat, DoF, f, n=None, m=None, xaxis = [-10,10,100], noplot=False):

    '''
    Plotting amplitude vs tuning for one LSC DoF.

    DoF    - Degree of freedom to sweep.
    f      - Frequency component relative to default
    n, m   - Mode numbers
    xaxis  - range to plot over [min, max, steps]
    '''
    
    _kat = kat.deepcopy()
    if isinstance(DoF, six.string_types):
        DoF = _kat.IFO.DOFs[DoF]
    # Adding detectors
    code = ""
    names = []
    for o in _kat.IFO.CAV_POWs:
        code += "{}\n".format(o.get_amplitude_cmds(f,n,m)[0])
        names.append(o.get_amplitude_name(f,n,m))
        
    # Adding simulation instructions
    code += pykat.ifo.scan_DOF_cmds(DoF, xlimits=[xaxis[0], xaxis[1]], steps=xaxis[2], relative=True)
    
    _kat.parse(code)
    out = _kat.run()
    
    if noplot:
        rtn = {'x': out.x}
        for n in names:
            rtn[n] = out[n]
        return rtn
    else:
        FS = 13
        LS = 12
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for n in names:
            ax.semilogy(out.x, out[n], label = n)
            
        ax.set_ylabel('$\mathrm{Amplitude}\ [\sqrt{\mathrm{W}}]$', fontsize=FS)
        ax.set_xlabel('$\mathrm{{ {}\ tuning\ [deg] }} $'.format(DoF.name), fontsize=FS)
        ax.set_xlim(out.x.min(), out.x.max())
        ax.grid()
        ax.legend(loc=3, fontsize=LS)
        plt.show(fig)
        
        return fig, ax
        
def amps_vs_dofs(kat, f, n=None, m=None, xaxis = [-1,1,100]):
    '''
    Plotting amplitude vs tuning for all LSC DoFs.

    f      - Frequency component relative to default
    n, m   - Mode numbers
    xaxis  - range to plot over [min, max, steps]
    '''

    dic = {}
    for d in kat.IFO.LSC_DOFs:
        xax = [d.scale*xaxis[0], d.scale*xaxis[1], xaxis[2]]
        dic[d.name] = amps_vs_dof(kat, d, f, n=n, m=m, xaxis = xax, noplot=True)

    N = len(dic)
        
    FS = 13
    LS = 11
    
    fig = plt.figure(figsize=(17,8))
    axs = []
    for k in range(N):
        axs.append(fig.add_subplot(2,3,k+1))
    for ax, v in zip(axs,dic.keys()):
        for n in dic[v].keys():
            if n != 'x':
                ax.semilogy(dic[v]['x'], dic[v][n], label = n)
        ax.set_ylabel('$\mathrm{Amplitude}\ [\sqrt{\mathrm{W}}]$', fontsize=FS)
        ax.set_xlabel('$\mathrm{{ {}\ tuning\ [deg] }} $'.format(v), fontsize=FS)
        ax.set_xlim(dic[v]['x'].min(), dic[v]['x'].max())
        ax.grid()
        ax.legend(loc=3, fontsize=LS)
    plt.show(fig)



def pows_vs_dof(kat, DoF, xaxis = [-10,10,100], noplot=False):

    '''
    Plotting amplitude vs tuning for one LSC DoF.

    DoF    - Degree of freedom to sweep.
    xaxis  - range to plot over [min, max, steps]
    '''
    _kat = kat.deepcopy()
    if isinstance(DoF, six.string_types):
        DoF = _kat.IFO.DOFs[DoF]
        
    # Adding detectors
    code = ""
    names = []
    for o in _kat.IFO.CAV_POWs:
        code += "{}\n".format(o.get_signal_cmds()[0])
        names.append(o.get_signal_name())
        
    # Adding simulation instructions
    code += pykat.ifo.scan_DOF_cmds(DoF, xlimits=[xaxis[0], xaxis[1]], steps=xaxis[2], relative=True)
    
    _kat.parse(code)
    out = _kat.run()
    
    if noplot:
        rtn = {'x': out.x}
        for n in names:
            rtn[n] = out[n]
        return rtn
    else:
        FS = 13
        LS = 12
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for n in names:
            ax.semilogy(out.x, out[n], label = n)
            
        ax.set_ylabel('$\mathrm{Power\ [W]}$', fontsize=FS)
        ax.set_xlabel('$\mathrm{{ {}\ tuning\ [deg] }} $'.format(DoF.name), fontsize=FS)
        ax.set_xlim(out.x.min(), out.x.max())
        ax.grid()
        ax.legend(loc=3, fontsize=LS)
        plt.show(fig)
        
        return fig, ax

def pows_vs_dofs(kat, xaxis = [-1,1,100]):
    '''
    Plotting amplitude vs tuning for all LSC DoFs.

    xaxis  - range to plot over [min, max, steps]
    '''

    dic = {}
    for d in kat.IFO.LSC_DOFs:
        xax = [d.scale*xaxis[0], d.scale*xaxis[1], xaxis[2]]
        dic[d.name] = pows_vs_dof(kat, d, xaxis = xax, noplot=True)

    N = len(dic)
        
    FS = 13
    LS = 11
    
    fig = plt.figure(figsize=(17,8))
    axs = []
    for k in range(N):
        axs.append(fig.add_subplot(2,3,k+1))
    for ax, v in zip(axs,dic.keys()):
        for n in dic[v].keys():
            if n != 'x':
                ax.semilogy(dic[v]['x'], dic[v][n], label = n)
        ax.set_ylabel('$\mathrm{Power\ [W]}$', fontsize=FS)
        ax.set_xlabel('$\mathrm{{ {}\ tuning\ [deg] }} $'.format(v), fontsize=FS)
        ax.set_xlim(dic[v]['x'].min(), dic[v]['x'].max())
        ax.grid()
        ax.legend(loc=3, fontsize=LS)
    plt.show(fig)

    return fig, axs

def pow_lsc(kat,xaxis = [-1,1,100]):
    '''
    Plotting cavity powers and error signals in the same figures. Can be used
    to see if error signal zero crossings coincide with power peaks.
    '''
    pass
    
