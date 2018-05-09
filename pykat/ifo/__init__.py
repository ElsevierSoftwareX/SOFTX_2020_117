from __future__ import print_function

import pykat
import pykat.exceptions as pkex
from pykat import isContainer
import numpy as np
import inspect
import math
import six
from copy import deepcopy
from pandas import DataFrame
from scipy.optimize import brute
from scipy.optimize import fmin
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from scipy.misc import comb



class SensingMatrix(DataFrame):
    @property
    def _constructor(self):
        return SensingMatrix
        
    def display(self):
        df = self
        from tabulate import tabulate
        keys = list(df.keys())
        keys.insert(0, "")
        
        print(tabulate(df.apply(np.abs), headers=keys, floatfmt=".3g"))
        print()
        
        keys[0] = "[deg]"
        print(tabulate(df.apply(lambda x: np.angle(x,True)), headers=keys, floatfmt=".3g"))
    
    def phase_reference(self, DOF):
        ref = np.angle(self.loc[DOF])
        return self.apply(lambda row: row * np.exp(-1j * ref), axis=1)
    
    def radar_plot(self, detector, _ax=None):
        import matplotlib.pyplot as plt
        import re 
        
        df = self
        
        A = df[detector]
        
        if _ax is None:
            ax = plt.subplot(111, projection='polar')
        else:
            ax = _ax
            
        ax.set_theta_zero_location('E')

        r_lim = (np.log10(np.abs(A)).min()-1, np.log10(np.abs(A)).max())

        for _ in A.keys():
            theta = np.angle(A[_])
            r = np.log10(np.abs(A[_]))

            #ax.plot((theta,theta), (r_lim[0], r), label=re.sub("[\(\[].*?[\)\]]", "", _).strip(), lw=2)
            ax.plot((theta,theta), (r_lim[0], r), lw=2, label=_)

        ax.set_title(detector)
        ax.set_ylim(r_lim[0], r_lim[1])
        ax.legend(bbox_to_anchor=(0.5, -0.1), loc="upper center", ncol=3) #len(A.keys()))
        ax.set_rticks(np.arange(*np.round(r_lim)))
        ax.grid(True, alpha=0.5, zorder=-10)

        if _ax is None:
            plt.tight_layout()
            plt.show()
        
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
                
    if len(components) != 0:
        raise pkex.BasePyKatException("Cannot find component {}".format(components))
        
    return kat
    
def round_to_n(x, n):
    if not x: return 0
    power = -int(math.floor(math.log10(abs(x)))) + (n - 1)
    factor = (10 ** power)
    return round(x * factor) / factor
    
def vprint(verbose, *printstr, **kwargs):
    if verbose:
        print(*printstr, **kwargs)
        
def BS_optical_path(thickness, n=1.44963098985906, angle=45.0):
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
         
def make_list_copy(_l):
    """
    Utility function, takes a list of strings or single string
    and returns a copy of a list, e.g.
    "string" copy to ["string"]
    ["string1", "string2"] copy to ["string1", "string2"]
    """
    if not isinstance(_l, (list, tuple, np.ndarray)):
        _l = [_l]
        
        
    return _l[:] # copy the list, just to be save

def scan_to_precision(DOF, target_precision, minmax="max", phi=0.0, precision=90.0, debug=None, extra_cmds=None):
    """
    Scans a DOF for a kat object to maximise or minimise the DOF's signal. 
    
    DOF - DOF object of the kat object
    target_precision - how accurate to max/min to
    minmax - "max" or "min"
    phi - initial starting point
    precision - look in phi-precision to phi+precision for min/max
    debug - output extra information
    extra_cmds - Extra commands to include in simulation run
    
    Returns phi of min/max and precision reached
    """
    while precision > target_precision * DOF.scale:
        out = scan_DOF(DOF.kat, DOF, xlimits = [phi-1.5*precision, phi+1.5*precision], extra_cmds=extra_cmds)
        
        phi, precision = find_peak(out, DOF.port.name, minmax=minmax, debug=debug)
         
    return phi, precision


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
    import pykat.external.peakdetect as peak
    
    stepsize = out.x[1]-out.x[0]
    if debug:
        print("  stepsize (precision) of scan: {0:g}".format(stepsize))

    _max, _min = peak.peakdetect( out[detector],out.x, 1)
    
    if debug:
        import matplotlib.pyplot as plt
        
        plt.figure()
        plt.semilogy(out.x,out[detector])
        plt.show()
        
        print("max: ")
        print(_max)
        print("min: ")
        print(_min)
        
    if len(_max) == 0 and minmax == "max":
        raise pkex.BasePyKatException("No maximum peaks found in {}".format(detector))
        
    if len(_min) == 0 and minmax == "min":
        raise pkex.BasePyKatException("No minimum peaks found in {}".format(detector))
        
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


def scan_optics_string(_optics, _factors, _varName='scan', target="phi", linlog="lin", xlimits=[-100, 100], steps=200, axis=1,relative=False):
    optics=make_list_copy(_optics)
    factors=make_list_copy(_factors)
    
    if len(optics) != len(factors):
        raise pkex.BasePyKatException("you must provide a factor for each optics")

    if linlog not in ["lin", "log"]: 
        raise pkex.BasePyKatException("linlog must be 'lin' or 'log'")
        
    _tuneStr  = "var {} 0\n".format(_varName)
    _tuneStr += "set {0}re {0} re\n".format(_varName)
    
    if axis==1:
        _tuneStr += "xaxis {} {} {} {} {} {}\n".format(_varName, 're', linlog, xlimits[0], xlimits[1], steps)
    elif (axis==2 or axis==3): 
        _tuneStr += "x{}axis {} {} {} {} {} {}\n".format(axis, _varName, 're', linlog, xlimits[0], xlimits[1], steps)
    else:
        raise pkex.BasePyKatException("axis must be 1, 2 or 3")

    _putStr = ""
    
    for idx, o in enumerate(optics):

        _putStr += "func sc{0} = ({1}) * ({2}) * ${3}re\n".format(o,np.abs(factors[idx]),np.sign(factors[idx]),_varName)
        _putStr += "noplot sc{}\n".format(o)
            
        if (relative):
            _putCmd = "put*"
        else:
            _putCmd = "put"
            
        _putStr += "{0} {1} {2} $sc{1}\n".format(_putCmd, o, target)
        # _putStr = "\n".join([_putStr, "{} {} {} {}{}".format(_putCmd, o, target, _xStr, axis)])            
        
    _tuneStr += _putStr

    #print(_tuneStr)
    #print()
    return _tuneStr

def scan_f_cmds(DOF, linlog="log", lower=10, upper=5000, steps=100):
    name = "_%s" % DOF.name
    cmds = DOF.fsig(name, 1)
    
    cmds += "\nxaxis {0} f {1} {2} {3} {4}\n".format(name, linlog, lower, upper, steps)
    
    return cmds

def scan_f(kat, DOF, linlog="lin", lower=10, upper=5000, steps=100, verbose=False):
    kat = kat.deepcopy()
    kat.parse(scan_f_cmds(DOF, linlog, lower, upper, steps))
    kat.verbose = verbose
    
    if DOF.port is not None:
        kat.parse(DOF.transfer())
        
    return kat.run(cmd_args=["-cr=on"])
    
def scan_DOF_cmds(DOF, xlimits=[-100, 100], steps=200, relative=False):    
    return scan_optics_string(DOF.optics, 
                              DOF.factors,
                              "scan",
                              target=DOF._mirror_target(),
                              linlog="lin",
                              xlimits=xlimits,
                              steps=steps,
                              axis=1,
                              relative=relative)

def scan_optic_cmds(self, optics, factors, xlimits=[-100, 100], steps=200,relative=False):
    return scan_optics_string(optics,
                              factors,
                              "scan",
                              linlog="lin",
                              xlimits=xlimits,
                              steps=steps,
                              axis=1,
                              relative=relative)
                              
def scan_DOF(kat, DOF, xlimits=[-100, 100], steps=200, relative=False, extra_cmds=None): 
    kat = kat.deepcopy()
    kat.parse(scan_DOF_cmds(DOF, xlimits=xlimits, steps=steps, relative=relative))
    
    if DOF.port is not None:
        kat.parse(DOF.signal())
    
    if extra_cmds:
        kat.parse(extra_cmds)
        
    return kat.run(cmd_args=["-cr=on"])

def scan_optics(kat, _optics, _factors, target="phi", xlimits=[-100, 100], steps=200,relative=False, extra_cmds=None): 
    """
    Scans one or more optics (by scanning the `target` parameter).
    
    Parameters:
    optics: list of names of components to be tuned
    factors: list of scaling factors for the tuning for each optics, first element must be 1.0
    xlimits: limits of the scan
    steps: number of steps to use in scan
    
    Usage:
        scan_optics(kat, "PRM", 1)
        scan_optics(kat, ["ETMX", "ETMY", [1, -1])
    """
    kat = kat.deepcopy()

    optics=make_list_copy(_optics)
    factors=make_list_copy(_factors)

    kat.parse(scan_optic_cmds(_optics, _factors, target=target, xlimits=xlimits, steps=steps, relative=relative))
    
    if extra_cmds:
        kat.parse(extra_cmds)
        
    return kat.run(cmd_args=["-cr=on"])

def optical_gain(DOF_sig, DOF_det, f=1, useDiff=True, deriv_h=1.0e-8):
    """
    Returns W/rad for length sensing. Will throw an exception if it can't convert the
    units into W/rad.
    """

    kat = DOF_sig.kat.deepcopy()
    kat.removeBlock('locks', False)

    if useDiff:
        _sigStr = DOF_sig.dcsig(deriv_h)
        _detStr = DOF_det.signal()
        _detName = DOF_det.signal_name()
    else:
        _sigStr = DOF_sig.fsig("sig1", fsig=f)
        _detStr  = DOF_det.transfer()
        _detName = DOF_det.transfer_name()
    
    kat.parse(_sigStr)
    kat.parse(_detStr)
    
    kat.noxaxis = True
    
    if useDiff:
        kat.parse("yaxis lin abs")
    else:
        kat.parse("yaxis lin abs:deg")
    
    out = kat.run()

    if useDiff:
        return float(out[_detName])*180/np.pi # W/rad
    else:
        if DOF_sig.sigtype == "phase":
            return float(np.real(out[_detName])) # W/rad
        elif DOF_sig.sigtype == "z":
            k = 2*np.pi/kat.lambda0
            return float(np.real(out[_detName])) / k # W/(m*k) = W/rad
        else:
            raise pkex.BasePyKatException("Not handling requested sigtype for unit conversion")

def diff_DOF(DOF, target, deriv_h=1e-12, scaling=1):
    """
    Returns commands to differentiate with respect to the DOF motion.
    
    This is typically used to find the slope of error signals at DC.
    """
    
    if target == "z":
        # As we can't target the z position of the component directly
        # we need to target phi and scale things
        target = 'phi'
        
    rtn = ("var x 0\n"
           "diff x re\n"
           "deriv_h {deriv_h}\n"
           "set _dx x re\n").format(deriv_h=deriv_h)
    
    for o,f in zip(DOF.optics, DOF.factors):
        rtn += "func sc{0} = ({1}) * ({2}) * ({3}) * $_dx\n".format(o,np.abs(f),np.sign(f),scaling)
        rtn += "noplot sc{}\n".format(o)
        rtn += "put* {0} {1} $sc{0}\n".format(o,target)
        
    return rtn
    
def scan_demod_phase_cmds(pd_detectors, demod_phase=1, relative=False, xaxis=1, steps=100, xlimits=(-180, 180)):
    """
    For a given list of detectors this will return the commands
    to scan the demod phase.
    
    pd_detectors: list of photodiode detector names
    demod_phase: which demodulation phase should be tuned, e.g. 2 would tune phase2 parameter
    relative: If true, put* is used
    xaxis: 1 for xaxis, 2 for x2axis scan
    steps: Number of steps to scan over
    xlimits: Range of scan in deg
    """
    pd_detectors = make_list_copy(pd_detectors)
    
    if xaxis not in [1, 2]:
        raise pkex.BasePyKatException("xaxis value must be 1 or 2")
        
    if xaxis == 1: xaxis = ""
        
    rtn = ("var scan 0\n"
           "x%saxis scan re lin %g %g %i\n" % (str(xaxis), xlimits[0], xlimits[1], steps))

    if relative:
        cmd = "put*"
    else:
        cmd = "put"
        
    for _ in pd_detectors:
        rtn += "%s %s phase1 $x%i\n" % (cmd, _, demod_phase)
        
    return rtn
    
def optimise_demod_phase(_kat, DOF, ports, minimise=False, debug=False):
    """
    This will optimise the demodulation phase at each port
    provide the largest slope in the detector outputs with respect
    to the DOF.
    
    Returns list of optimised demodulation phases corresponding to the order
    of detectors given in pd_detectors
    """
    kat = _kat.deepcopy()
    
    if isinstance(DOF, six.string_types):
        DOF = kat.IFO.DOFs[DOF]
        
    # Get a list of port objects even if user provided them in string names
    _ports = []
    for port in ports:
        if isinstance(port, six.string_types):
            _ports.append(kat.IFO.Outputs[port])
        else:
            _ports.append(kat.IFO.Outputs[port.name])
            
        if _ports[-1].f is None:
            raise pkex.BasePyKatException("port %s cannot have its demodulation phase optimised as it isn't demodulated" % port.name)
    
    ports = _ports
    
    pd_detectors = []
        
    # Add in the signals
    for port in ports:
        pd_detectors.append(kat.IFO.Outputs[port.name].add_signal("I", DOF.sigtype))
    
    if debug:
        print("Optimising pds: %s" % pd_detectors)
        
    kat.removeBlock("locks", False)
    kat.removeBlock("powers", False)
    
    kat.parse( diff_DOF(DOF, DOF.sigtype), addToBlock="OPTIMISE")
    kat.parse( scan_demod_phase_cmds(pd_detectors), addToBlock="OPTIMISE")
    
    # Analyitcally we can find the phase which gives a maxmium
    # by solving for A and B in:
    #   f(x) = Acos(x+B)
    # assuming we take the two data points x = {0, pi/2} we find
    #   B = arctan(y2/y1)
    kat.xaxis.limits = (0, 90)
    kat.xaxis.steps = 1
    
    if debug:
        print(kat & "OPTIMISE")
        
    out = kat.run(cmd_args=["-cr=on"])
    
    if debug:
        print(out.y)
        
    rtn = []
    
    for pd, port in zip(pd_detectors, ports):
        
        y1, y2 = out[pd]
        x = np.deg2rad(out.x)
        R = np.sqrt(y1**2 + y2**2)
        phi = np.rad2deg(np.arctan2(y2,y1))
        
        if minimise:
            phi += 90
        
        # All in I quadrature so no
        _kat.IFO.Outputs[port.name].phase = phi
        
        rtn.append(phi)
        
    
    return rtn
    
def mismatch_cavities(base, node):
    _kat = base.deepcopy()
    
    _kat.removeBlock("locks", False)
    _kat.removeBlock("errsigs", False)

    for _ in list(_kat.detectors.values()):
        _.remove()

    _kat.parse("bp qx x q "+node)
    _kat.parse("bp qy y q "+node)
    _kat.noxaxis = True
    _kat.yaxis = "re:im"
    _kat.maxtem = 0
    
    qs  = []
    qxs = []
    qys = []

    gauss = list(_kat.getAll(pykat.commands.cavity))
    
    for _ in _kat.nodes.getNodes():
        if _kat.nodes[_].q is not None:
            gauss.append(_kat.nodes[_])
    
    for q in gauss:
        # Switch off all cavities
        for _ in gauss:
            _.enabled = False
            
        # Then select one at a time
        q.enabled = True
        out = _kat.run()

        qs.append(q.name)
        
        qxs.append(pykat.BeamParam(q=out['qx']))
        qys.append(pykat.BeamParam(q=out['qy']))

    mmx = DataFrame(index=qs, columns=qs)
    mmy = DataFrame(index=qs, columns=qs)
    
    for c1, q1x, q1y in zip(qs, qxs, qys):
        for c2, q2x, q2y in zip(qs, qxs, qys):
            mmx[c1][c2] = pykat.BeamParam.overlap(q1x, q2x)
            mmy[c1][c2] = pykat.BeamParam.overlap(q1y, q2y)
            
    return 1-mmx.astype(float), 1-mmy.astype(float), list(zip(qs, qxs, qys))

def mismatch_scan_RoC(base, node, mirror, lower, upper, steps):
    _kat = base.deepcopy()
    _kat.removeBlock("locks", False)
    _kat.removeBlock("errsigs", False)

    for _ in list(_kat.detectors.values()):
        _.remove()

    _kat.parse("bp qx x q "+node)
    _kat.parse("bp qy y q "+node)
    _kat.yaxis = "re:im"
    _kat.maxtem = 0
    _kat.parse("""
    xaxis* {name} Rcx lin -1 1 10
    put {name} Rcy $x1
    """.format(name=mirror))
    
    _kat.xaxis.limits = (lower, upper)
    _kat.xaxis.steps = steps

    out = _kat.run()
    return out['qx'], out['qy']
    
def mismatch_scan_RoC_2D(base, node, mirror1, lower1, upper1, steps1, mirror2, lower2, upper2, steps2):
    _kat = base.deepcopy()
    _kat.removeBlock("locks", False)
    _kat.removeBlock("errsigs", False)

    for _ in list(_kat.detectors.values()):
        _.remove()

    _kat.parse("bp qx x q "+node)
    _kat.parse("bp qy y q "+node)
    _kat.yaxis = "re:im"
    _kat.maxtem = 0 # don't need any maxtem as we are just tracing
    _kat.parse("""
    xaxis* {name1} Rcx lin 0 1 10
    x2axis* {name2} Rcx lin 0 1 10
    put {name1} Rcy $x1
    put {name2} Rcy $x2
    """.format(name1=mirror1, name2=mirror2))
    
    _kat.xaxis.limits = (lower1, upper1)
    _kat.xaxis.steps = steps1
    
    _kat.x2axis.limits = (lower2, upper2)
    _kat.x2axis.steps = steps2
    _kat.retrace = "force"
    
    out = _kat.run()
    
    return out['qx'], out['qy']

def mismatch_scan_L(base, node, length, lower, upper, steps):
    _kat = base.deepcopy()
    _kat.removeBlock("locks", False)
    _kat.removeBlock("errsigs", False)

    for _ in list(_kat.detectors.values()):
        _.remove()

    _kat.parse("bp qx x q "+node)
    _kat.parse("bp qy y q "+node)
    _kat.yaxis = "re:im"
    _kat.maxtem = 0
    _kat.parse("""
    xaxis* {name} L lin -1 1 10
    """.format(name=length))
    
    _kat.xaxis.limits = (lower, upper)
    _kat.xaxis.steps = steps

    out = _kat.run()
    return out['qx'], out['qy']


def modematch(kat, components, cavs, node, verbose = False):
    '''
    Mode matches the cavity eigenmmodes for the cavities in cavs by varying the
    components in components. Computes mode overlaps between the cavity eigenmodes.
    Minimises the maximum cavity mismatch. 
    
    Inputs
    --------
    kat         - kat-object to run
    components  - list of names of components to vary. Each entry must be a lens (varies f), 
                  mirror (Rc), BS (Rc) or space (L). 
    cavs        - list with names of cavities to match
    node        - name of node where to compute mismatches. Must be first or
                  last node in IFO for reliable result.
    verbose     - If true, prints the new optimised paramaters
    
    Returns
    --------
    kat1        - kat-object with the optimised component parameters. Deepcopy of kat.
    out         - array with the new optimised values
    
    '''
    Nc = len(cavs)
    Np = len(components)
    kat1 = kat.deepcopy()
    kat2 = kat.deepcopy()
    # Storing parameters to tune, and their initial values, in lists
    attrs = []
    attrs2 = []
    p0 = []
    for c in components:
        if isinstance(kat1.components[c], pykat.components.lens):
            p0.append(kat1.components[c].f.value)
        elif (isinstance(kat1.components[c], pykat.components.mirror) or 
              isinstance(kat1.components[c], pykat.components.beamSplitter)):
            p0.append(kat1.components[c].Rc.value)
        elif isinstance(kat1.components[c], pykat.components.space):
            p0.append(kat1.components[c].L.value)

        attrs.append(kat1.components[c])
        attrs2.append(kat2.components[c])

            
    # Switching off cavity commands for cavities not in cavs
    for cav in kat1.getAll(pykat.commands.cavity):
        if not cav.name in cavs:
            cav.remove()

    # Cost function
    def func(p):
        for k in range(Np):
            if isinstance(attrs[k], pykat.components.lens):
                attrs[k].f = p[k]
            elif isinstance(attrs[k], pykat.components.space):
                attrs[k].L = p[k]
            elif (isinstance(attrs[k], pykat.components.mirror) or 
                  isinstance(attrs[k], pykat.components.beamSplitter)):
                attrs[k].Rc = p[k]
        
        mmx, mmy, qs = pykat.ifo.mismatch_cavities(kat1, node)
        mm = np.zeros(comb(Nc,2,exact=True), dtype=float)
        cs = deepcopy(cavs)
        k = 0
        for c1 in cavs:
            cs.pop(0)
            for c2 in cs:
                mm[k] = np.sqrt(mmx[c1][c2])*np.sqrt(mmy[c1][c2])
                k += 1
        # print(kat1.CPN_TL.f.value, kat1.CPW_TL.f.value, mm.mean())
        return mm.max()
        
    if verbose: 
        # Computing initial mismatch. Only for display
        mmx, mmy, qs = pykat.ifo.mismatch_cavities(kat1, node)
        mm0 = np.zeros(comb(Nc,2,exact=True), dtype=float)
        cs = deepcopy(cavs)
        k = 0
        for c1 in cavs:
            cs.pop(0)
            for c2 in cs:
                mm0[k] = np.sqrt(mmx[c1][c2])*np.sqrt(mmy[c1][c2])
                k += 1
        
    # Optimising
    opts = {'xtol': 1.0, 'ftol': 1.0e-7, 'disp': False}
    out = minimize(func, p0, method='Nelder-Mead', options=opts)
    
    if not out['success']:
        pkex.printWarning(out.message)
        
    # Setting new parameters to kat-object
    for k in range(Np):
        if isinstance(attrs2[k], pykat.components.lens):
            attrs2[k].f = out.x[k]
        elif isinstance(attrs2[k], pykat.components.space):
            attrs2[k].L = out.x[k]
        elif (isinstance(attrs2[k], pykat.components.mirror) or 
              isinstance(attrs2[k], pykat.components.beamSplitter)):
            attrs2[k].Rc = out.x[k]

    if verbose:
        print('Maximum mismatch: {:.2e} --> {:.2e}'.format(mm0.max(), out.fun))
        for c in components:
            if isinstance(kat.components[c], pykat.components.lens):
                print(' {}.f: {:.5e} m --> {:.5e} m'.format(c, kat.components[c].f.value, kat2.components[c].f.value ))
            elif isinstance(kat.components[c], pykat.components.space):
                print(' {}.L: {:.5e} m --> {:.5e} m'.format(c, kat.components[c].L.value, kat2.components[c].L.value ))
            elif (isinstance(kat.components[c], pykat.components.mirror) or 
                  isinstance(kat.components[c], pykat.components.beamSplitter)):
                print(' {}.Rc = {:.5e} m --> {:.5e} m'.format(c, kat.components[c].Rc.value, kat2.components[c].Rc.value ))
                      
    return kat2, out.x



def c2r_errsig(z, phase):
    '''
    Function that returns the real error signal by setting a demodulation phase
    to the complex error signal.

    Paramaters:
    -----------
    z      - array containing the complex error signal created by I + 1j*Q, where I and Q are
             the real error signals with with 0 deg and 90 deg demodulation phases, respectively. 
    phase  - the demodulation phase [deg].
    '''
    return np.real(z*np.exp(-1j*phase*np.pi/180.0))

def opt_demod_phase(cdata, x, xbounds=None, err_tol=1e-5, xatol=1e-9, isplot=False):
    '''
    Optimizes the demodulation phase of a complex error signal generated by Finesse.
    Demands that the error signal is smaller than err_tol to compute the slope or
    optical gain. 

    Paramaters:
    -----------
    cdata    - Array with a complex error signal from Finesse.
    x        - Array with DoF values. Typically kat.run().x from Finesse. 
    xbounds  - x-range where to search.
    err_tol  - Defines how small the error signal must be to be considered
               a zero-crossing
    xatol    - accuracy in x
    isplot   - Plot error signal or not.

    Returns:
    -----------
    demod_phase   - The optimal demodulation phase [deg]
    optical_gain  - The optical gain (slope of the error signal). The unit
                    depends on the unit of x.
    '''

    if xbounds is None:
        xbounds = np.array([x[0],x[-1]])
    xscale = xbounds[-1]-xbounds[0]
    xn = x/xscale
    xbounds = xbounds/xscale
    yscale = 2*np.abs(cdata).max()
    cdatan = cdata/yscale
    
    # Function to optimise
    def max_og(phi):
        err = c2r_errsig(cdatan, phi)
        dErr = np.gradient(err, xn[1]-xn[0])
        f = interp1d(xn, err, kind='linear')
        df = interp1d(xn, dErr, kind='quadratic')

        def find_op(x1):
            y1 = np.abs(f(x1))
            dy1 = df(x1)
            # print(phi, y1, dy1)
            if y1 > err_tol:
                if dy1 <= 0:
                    return np.inf
                return y1
            else:
                #print(dy1, phi)
                # print(x,y1,dy1)
                #print('hej!')
                return -dy1
        sol_2 = minimize_scalar(find_op, bounds=xbounds, method='bounded', 
                                options={'maxiter': 500, 'disp': True, 'xatol': xatol})
        # print(sol_2)
        return sol_2.fun

    # Demodulation phases to look at with brute force. 
    min_phase = -180
    max_phase = 180
    step = 1
    demod_phase_bounds = slice(min_phase, max_phase, step)
    phase_bounds = (demod_phase_bounds,)
    # Searching for optimal demodulation phase.
    sol = brute(max_og, ranges=phase_bounds, full_output=True, finish=fmin, disp=False)
    # Reading out results
    demod_phase = sol[0][0]
    optical_gain = -sol[1]*yscale/xscale


    if optical_gain<=0:
        print(('Error: Found no error signal small enough in the given domain. '
                   'Change the domain, and/or decrease step sizes'))
        
        # PLOTTING WILL BE REMOVED
        if isplot:
            import matplotlib.pyplot as plt
            phases = np.linspace(min_phase, max_phase, (max_phase - min_phase)/step + 1)
            cm = plt.cm.Spectral_r
            norm = mpl.colors.Normalize(vmin=phases[0], vmax=phases[-1])
            s_m = mpl.cm.ScalarMappable(cmap=cm, norm=norm)
            s_m.set_array([])
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for p in phases[::8]:
                c = s_m.to_rgba(p)
                ax.plot(x,c2r_errsig(cdata,p),color=c)
            ax.set_xlim(x.min(),x.max())
            cb = plt.colorbar(s_m)

            cb.set_label('Demod. phase [deg]')
            plt.show(fig)
        return None, None

    return demod_phase, optical_gain

def ASC_demod_phase(kat, asc_dof, output, xaxis=[-1e-6, 1e-6, 50], err_tol = 1e-5,
                    xatol = 1e-9, isplot=False, verbose=False):

    _kat = kat.deepcopy()
    cmd_str = output.get_signal_cmds(quad='I',sigtype=asc_dof.sigtype)[0].split()
    cmd_str[1] = cmd_str[1][:-4] + cmd_str[1][-2:]
    _kat.parse("""
    {0[0]} {0[1]} {0[2]} {0[4]}
    pdtype {0[1]} y-split
    """.format(cmd_str), addToBlock=output._block)

    _kat.parse(pykat.ifo.scan_DOF_cmds(asc_dof,xlimits=[xaxis[0],xaxis[1]],steps=xaxis[2],relative=False))
    _kat.parse('yaxis abs:deg')

    out = _kat.run()
    y = out[cmd_str[1]]
    x = out.x

    demod_phase, og = opt_demod_phase(y, x, xbounds=xaxis[:-1], err_tol=err_tol, xatol=xatol, isplot=isplot)

    # THE PLOTTING OPTION WILL BE REMOVED
    if isplot:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x, c2r_errsig(y,demod_phase), 'b-', label='{} at {}'.format(asc_dof.name,cmd_str[1]))
        ax.set_xlabel("{} {} [rad]".format(asc_dof.name, asc_dof.sigtype))
        ax.set_ylabel("{} at {} [W]".format(asc_dof.name, cmd_str[1]))
        ax.set_xlim(x.min(), x.max())
        ax.legend(loc=2)
        ax.grid()
        plt.show(fig)
    if verbose:
        print('{} at {}: demod. phase =  {:.3e} deg, optical gain = {:.3e} W/rad'
              .format(asc_dof.name,cmd_str[1],demod_phase, og))

    return demod_phase, og
    

def LSC_demod_phase(lsc_dof, xaxis=[-.1, .1, 50], err_tol = 1e-5, pwr_dof=None,
                    P_tol = 1e-5, xatol = 1e-9, isMax = True, isplot=False, verbose=False):
    '''
    Finding the demodulation phase of a LSC DOF that generates the maximum slope of the error signal.

    Parameters:
    -----------
    kat      - Kat object to run
    lsc_dof  - pykat.ifo.DOF object, for which the error signal demodulation phase is optimised
    xaxis    - Relative range around the current operating point, in which we demand the zero
               crossing of the error signal to occur.
    err_tol  - The error singal is demanded to be smaller than err_tol at the point where the
               slope is computed.
    pwr_dof  - Optional power degree of freedom that defines a power in the relevant cavity. If set,
               the optimal tuning for the power is first found, and this is used to specify the
               relative range in which the error signal zero crossing must occur.
    P_tol    - How much the operating point is allowed to deviate from maximum/minimum power.
    xatol    - Absolute tolerance in tuning when finding the maximum slope.
    isMax    - If true, the operating point should be at a power peak, otherwise at a power minimum.
    isPlot   - Plot result.
    verobse  - Print result.

    Returns:
    -----------
    demod_phase   - The optimal demodulation phase [deg]
    og            - The optical gain (slope of the error signal [W/deg])
    '''
    # THE PLOTTING OPTION WILL BE REMOVED

    _kat = lsc_dof.kat.deepcopy()

    # Scan instructions
    scstr = pykat.ifo.scan_DOF_cmds(lsc_dof, xlimits=[xaxis[0], xaxis[1]], steps=xaxis[2], relative=True)
    # LSC detector
    lscstr = lsc_dof.signal()[0]        

    # If pwr_dof specified, finding power peak and a tuning range around it where the power
    # change is < P_tol. If not pwr_dof is specifed we set the range to the xaxis range. 
    if pwr_dof is not None:
        # Power detector
        pwrstr = pwr_dof.signal()[0]
        [Pmax, xmax], Pr, P_func = find_max_power(_kat, scstr+pwrstr, pwr_dof.signal_name(),
                                                  P_tol = P_tol, isplot=False, isMax=isMax)
    else:
        Pr = np.array(xaxis[0:2])

    # Checking if demodulation is used.
    if len(lscstr.split()) < 3:
        print('Error: No demodulation in detector {}'.fomrat(lsc_dof.signal_name()))
        return None

    # Changing detector and y-axis to obtain complex error signal.
    lscstr = lsc_dof.signal()[0].split()
    lscstr = "{0[0]} {0[1]} {0[2]} {0[4]}".format(lscstr)
    _kat.yaxis = 'abs:deg'

    # Parsing and setting simulation instructions.
    _kat.parse(scstr+lscstr)
    _kat.xaxis.limits = np.array(Pr)
    # Running
    out = _kat.run()
    # Extracting results
    y = out[lsc_dof.signal_name()]
    x = out.x

    demod_phase, og = opt_demod_phase(y, x, xbounds=Pr, err_tol=err_tol, xatol=xatol, isplot=isplot)

    # PLOTTING WILL BE REMOVED
    if isplot:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x, c2r_errsig(y,demod_phase), 'r-', label='{} error signal'.format(lsc_dof.name))
        ax.set_xlabel("{} tuning [deg]".format(lsc_dof.name))
        ax.set_ylabel("{} err. sig [W]".format(lsc_dof.name), color='r')
        ax.set_xlim(x.min(), x.max())
        ax.legend(loc=2)
        if pwr_dof is not None:
            ax2 = ax.twinx()
            ax2.plot(x, P_func(x), 'b-', label='Power')
            ax2.set_ylabel('{} [W]'.format(pwr_dof.port.name), color='b')
            ax2.set_xlim(x.min(), x.max())
            ax2.legend(loc=4)
        ax.grid()
        plt.show(fig)
    if verbose:
        print('{}: demod. phase =  {:.3e} deg, optical gain = {:.3e} W/deg'.format(lsc_dof.name,demod_phase, og))
    return demod_phase, og

def find_max_power(kat, scanstring, detector, P_tol=1e-4, isplot=False, isMax = True):
    '''
    Finds the maximum power in the detector, in the run defined by the
    kat object + scanstring
    '''
    _kat = kat.deepcopy()
    _kat.parse(scanstring)
    out = _kat.run()
    x0 = out.x
    P0 = out[detector]
    if isMax:
        i = P0.argmax()
    else:
        i = P0.argmin()

    #print(P0[i])
    #print(out[pwr_dof.signal_name()])

    if isplot:
        import matplotlib.pyplot as plt
        plt.plot(x0, P0)
        plt.xlim(x0.min(),x0.max())
        plt.grid()
        plt.show()

    # tune = pykat.ifo.find_peak(out, pwr_dof.signal_name(), minmax='max', debug=False)
    
    _kat.xaxis.limits = [out.x[i-1], out.x[i+1]]
    out = _kat.run()
    x = out.x
    P = out[detector]

    if isMax:
        f = interp1d(x, -P, kind='quadratic')
    else:
        f = interp1d(x, P, kind='quadratic')
        
    sol = minimize_scalar(f, method='bounded', bounds = (x0[i-1], x0[i+1]), tol=None,
                          options={'xatol': 1.48e-8, 'maxiter': 500})
    if isMax:
        Pmax = -sol.fun
    else:
        Pmax = sol.fun
    xmax = sol.x

    f2 = interp1d(x0, P0, kind='quadratic')
    def f3(x):
        if isMax:
            return np.abs(f2(x) - (1-P_tol)*Pmax)
        else:
            return np.abs(f2(x) - (1+P_tol)*Pmax)

    
    #sol = minimize_scalar(f, method='brent', tol=None,
    #                      options={'xtol': 1.48e-8, 'maxiter': 500})

    sol = minimize_scalar(f3, method='bounded', bounds = (x0[i-1], x0[-1]), tol=None,
                          options={'xatol': 1.48e-8, 'maxiter': 500})
    
    x1 = sol.x
    dx = np.abs(x1 - xmax)
    x1 = xmax - dx
    x2 = xmax + dx
    Prange = [x1, x2]
    
    # print(xmax)
    # print(dx,x1,x2)
    
    if isplot:
        plt.plot(x, P)
        if isMax:
            plt.plot([xmax, xmax], [P.min(), Pmax],'r--')
            plt.plot([x1,x1], [P.min(), Pmax],'g--')
            plt.plot([x2,x2], [P.min(), Pmax],'g--')
        else:
            plt.plot([xmax, xmax], [Pmax, P.max()], 'r--')
            plt.plot([x1,x1], [Pmax, P.max()], 'g--')
            plt.plot([x2,x2], [Pmax, P.max()],'g--')
            
        plt.plot([x.min(), x.max()], [Pmax, Pmax],'r--')
        plt.plot([x.min(), x.max()], [(1+P_tol)*Pmax, (1+P_tol)*Pmax],'g--')
        plt.xlim(x.min(),x.max())
        plt.grid()
        plt.show()
        

    return [Pmax, xmax], Prange, f2









class IFO(object):
    """
    A generic object that contains various interferometer properties.
    
    This should contain methods that operate on the IFO or underlying kat object
    in a detector agnostic manner. e.g. scanning a particular DOF.
    
    Functions that directly alter the kat object, such as setting lengths, should
    be contained in a detector specific class that inherits from this, such as
    ALIGO_IFO.
    """
    def __init__(self, kat, tuning_keys_list, tunings_components_list):
        self.__kat = kat
        self.__tuning_keys = frozenset(make_list_copy(tuning_keys_list))
        self.__tuning_comps = make_list_copy(tunings_components_list[:])
    
    @property
    def kat(self): return self.__kat
    
    def requires_detectors(self, *args):
        if len(args) == 1 and pykat.isContainer(args[0]):
            self.requires_detectors(*args[0])
        else:
            for _ in args:
                if _ not in self.kat.detectors:
                    raise pkex.BasePyKatException("Detector `%s` was not found in the kat object" % _)
    
    def requires_DOFs(self, *args):
        if len(args) == 1 and pykat.isContainer(args[0]):
            self.requires_DOFs(*args[0])
        else:
            for _ in args:
                if _ not in self.DOFs:
                    raise pkex.BasePyKatException("DOF `%s` was not found in the IFO object" % _)

    def _tuning_key(self, **kwargs):
        if set(kwargs.keys()) != self.__tuning_keys:
            raise pkex.BasePyKatException("input keyword arguments should be: %s" % ", ".join(self.__tuning_keys))
            
        vals = []
        
        for key in self.__tuning_keys:
            vals.append(str(kwargs[key]))
        
        return "-".join(vals)
        
    def zero_locks(self, verbose=False):
        """
        Zeroes the error signals by running the currently setup locks in this kat object.
        The corrected tunings are then applied to this object.  
        """
        
        if verbose:
            print("Old tunings")
            print(self.get_tunings())
        
        base = self.kat.deepcopy()
        base.noxaxis = True
        base.verbose = False
        out = base.run(cmd_args=["-cr=on"])
        self.apply_lock_feedback(out)
        
        if verbose:
            print("New tunings")
            print(self.get_tunings())
            
    def get_tuning_comps(self):
        return self.__tuning_comps

            
    def get_tunings(self):
        """
        For the current state of the kat object, this will return 
        a dictionary that contains all the required tunings for the interferometer.
        """
        keys = self.__tuning_keys
            
        rtn = {}
        
        for comp in self.__tuning_comps:
            if not hasattr(self.kat, comp):
                raise pkex.BasePyKatException("`%s` is not a component of the kat object" % comp)
                
            _ = getattr(self.kat, comp)
            
            rtn[comp] = float(_.phi)
        
        rtn["keys"] = {}
        
        for key in keys:
            rtn["keys"][key] = getattr(self.kat, key)
        
        return rtn
            
    def save_tunings(self):
        tunings = self.kat.IFO.get_tunings()
        
        if "IFO.tunings" not in self.kat.data:
            self.kat.data["IFO.tunings"] = {}
            
        key = self._tuning_key(**tunings['keys'])
        self.kat.data["IFO.tunings"][key] = tunings
        
        return tunings
    
    def load_tunings(self, **kwargs):
        """
        sets the tunings of optical components and the corresponding values
        for maxtem and phase
        """
        
        key = self._tuning_key(**kwargs)
        
        if "IFO.tunings" not in self.kat.data:
            pkex.printWarning("No IFO tunings are stored in this kat")
        elif key not in self.kat.data["IFO.tunings"]:
            pkex.printWarning("Could not find tunings for values %s" % kwargs)
        else:        
            return self.kat.data["IFO.tunings"][key]
    
    def apply_tunings(self, tunings):
        for comp in tunings:
            if comp == "keys": continue
            
            if comp in self.kat.components:
                self.kat.components[comp].phi = tunings[comp]
            else:
                pkex.printWarning("%s not present in kat, skipping" % comp)
                
    def strsToDOFs(self, DOFs):
        """
        Converts a list of strings of DOF names into a
        list of DOF objects. Throws errors if DOF name 
        are not found.
        """
        
        dofs = []
        
        for _ in DOFs:
            if isinstance(_, six.string_types):
                if _ in self.DOFs:
                    dofs.append(self.DOFs[_])
                else:
                    raise pkex.BasePyKatException("Could not find DOF called `%s`. Possible DOF options: %s" % (_, str(list(self.DOFs.keys()))))
            else:
                raise pkex.BasePyKatException("'%s' not possible DOF options: %s" % (_, str(list(self.DOFs.keys()))))
        
        return dofs
    
    def sensing_matrix(self, DOFs, detectors, frequency=1):
        """
        Computes a sensing matrix for a collection of DOFs and detectors at a particular signal frequency.
        
        The function returns an augmented Pandas DataFrame object. This function adds two methods to the
        DataFrame:
            sens = kat.IFO.sensing_matrix(DOFs, detectors, 1)
            sens.print() # Will show an ascii table of the sensing matrix
            sens.radar_plot(detector_name) # Will display sensing matrix in polar plot form
        
        DOFs: collection of DOF objects
        detectors: list of detector names
        frequency: frequency to compute sensing matrix at [Hz]
        Returns: Pandas DataFrame
        """
        self.requires_detectors(detectors)
        self.requires_DOFs(DOFs)
        
        if not isContainer(DOFs): DOFs = [DOFs]
        if not isContainer(detectors): detectors = [detectors]
        
        data = []

        for DOF in DOFs:
            kat = self.kat.deepcopy()
    
            kat.noxaxis = True
            kat.yaxis = "re:im"
            kat.removeBlock("locks", False)
            kat.removeBlock("powers", False)
            kat.removeBlock("errsigs", False)
    
            kat.parse(self.DOFs[DOF].fsig(fsig=frequency) )
    
            data.append(kat.run()[detectors])
    
        return SensingMatrix(data, columns=detectors, index=DOFs)
    
class DOF(object):
    """
    Defining a degree of freedom for the interferometer, includes the
    objects and how to move them, and the default output port to read
    out the DOF signal.
    """
    def __init__(self, IFO, _DOFName, _port, _quad, _optics, _factors, _scale, sigtype="z"):
        self.__IFO = IFO
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
    
    def _mirror_target(self):
        """
        Returns which parameter to target in puts and xaxis depending on sigtype
        """
    
        if self.sigtype == "z":
            target = "phi"
        elif self.sigtype == "pitch":
            target = "ybeta"
        elif self.sigtype == "yaw":
            target = "xbeta"
        else:
            raise pkex.BasePyKatException("Unexpected sigtype %s" % self.sigtype)
        
        return target
        
    @property
    def kat(self):
        """For referencing the kat object this DOF is associated with"""
        return self.__IFO.kat
        
    def scan(self, **kwargs):
        """
        Convenience method for calling `scan_DOF` for this particular DOF and associated kat object.
        
        See `pykat.ifo.scan_DOF` for keyword arguments options.
        """
        return scan_DOF(self.__IFO.kat, self, **kwargs)
        
    def scan_f(self, *args, **kwargs):
        """
        Runs an fsig simulation scaning this DOF
        
        See `pykat.ifo.scan_f` for keyword arguments options.
        """
        return scan_f(self.__IFO.kat, self, *args, **kwargs)
        
    def apply_tuning(self, phi, add=False):
        for idx, o in enumerate(self.optics):
            if add:
                self.__IFO.kat.components[o].phi += phi * self.factors[idx]
            else:
                self.__IFO.kat.components[o].phi = phi * self.factors[idx]
    
    def add_signal(self):
        if self.port is None:
            raise pkex.BasePyKatException("No port is associated with {}".format(self.name))
        
        return self.port.add_signal(self.quad, self.sigtype)
                
    def signal(self):
        if self.port is None:
            raise pkex.BasePyKatException("No port is associated with {}".format(self.name))
            
        return self.port.get_signal_cmds(quad=self.quad, sigtype=self.sigtype)
        
    def signal_name(self):
        if self.port is None:
            raise pkex.BasePyKatException("No port is associated with {}".format(self.name))
            
        return self.port.get_signal_name(self.quad, self.sigtype)

    def transfer(self, phase2=None):
        if self.port is None:
            raise pkex.BasePyKatException("No port is associated with {}".format(self.name))
            
        return self.port.get_transfer_cmds(self.quad, phase2=phase2)
        
    def transfer_name(self):
        if self.port is None:
            raise pkex.BasePyKatException("No port is associated with {}".format(self.name))
            
        return self.port.get_transfer_name(self.quad)

    def fsig(self, _fsigName=None, fsig=1.0):
        _fsigStr= ""
        
        if _fsigName is None:
            _fsigName = self.name+ "_fsig"
        
        for idx, o in enumerate(self.optics):
            phase = 0.0
            
            if self.factors[idx] < 0:
                phase = 180.0
            
            _fsigStr = "\n".join([_fsigStr, "fsig {name} {component} {type} {f} {phase} {amp}".format(name=_fsigName, 
                                                                                                    component=o,
                                                                                                    type=self.sigtype,
                                                                                                    f=fsig,
                                                                                                    phase=phase,
                                                                                                    amp=abs(self.factors[idx]))])
            
        return _fsigStr

    def dcsig(self, deriv_h=1e-8):
        '''
        Returns Finesse code for computing the DC slope of the error signal by
        using the Finesse command diff. 
        '''
        return diff_DOF(self, self._mirror_target(), deriv_h=deriv_h)

class Output(object):
    """
    This object defines a location in an interferometer where detectors are placed and demodulated at a particular
    modulation frequency or DC. It does not specify any detectors in particular.
    However, using this object you can add the following detectors to the associated kat object:

        * Photodiodes, for error signals (signal)
        * Transfer functions to use with fsig (transfer)
        * Amplitude detectors (amplitude)
    
    In brackets are the tags associated with each type of detector. The functions
    here are named with each tag. The possible options are:
    
    You can add many detectors at a given port, which readout different quadratures or types of transfer functions
    or signals.
    """
    def __init__(self, IFO, name, nodeNames, f_property_name=None, phase=0, block=None):
        if f_property_name is not None and not isinstance(f_property_name, six.string_types):
            raise Exception("f_property_name should be a property name of a frequency in class {}".format(IFO.__class__))
            
        self.__IFO = IFO
        self.name = name
        self.nodeNames = make_list_copy(nodeNames)
        self.__f_property_name = f_property_name
        self.phase = phase    # demodulation phase for I quadrature, float
        self._block = block
    
    @property
    def kat(self):
        """For referencing the kat object this DOF is associated with"""
        return self.__IFO.kat
          
    @property
    def f(self):
        if self.__f_property_name is None:
            return None
        
        return getattr(self.__IFO, self.__f_property_name)
        
    def check_nodeName(self):
        self.nodeName = None
        
        for node in self.nodeNames:
            _node = node
            if _node[-1] == "*":
                _node = _node[:-1]
                
            if _node in self.__IFO.kat.nodes:
                self.nodeName=node
                break
                
        if self.nodeName==None:
            raise pkex.BasePyKatException("port {}: cannot find any of these nodes: '{}'".format(self.name,self.nodeNames))

    def get_amplitude_name(self, f, n=None, m=None):
        ''''
        Returning name of amplitude detector
        '''
        if round(abs(f)) < 1e3:
            fstr = "{}".format(round(f))
        elif round(abs(f)) < 1e6:
            fstr = "{}k".format(round(f/1e3))
        elif round(abs(f)) < 1e9:
            fstr = "{}M".format(round(f/1e6))
        elif round(abs(f)) < 1e12:
            fstr = "{}G".format(round(f/1e9))
        elif round(abs(f)) < 1e15:
            fstr = "{}T".format(round(f/1e12))
        else:
            fstr = "{:.0e}".format(f)
            
        if n is None or m is None:
            rtn = "{}_{}_ad".format(self.name,fstr)
        else:
            rtn = "{}_{}_{}{}_ad".format(self.name,fstr,n,m)
        return rtn
    
    def get_amplitude_cmds(self, f, n=None, m=None):
        rtn = []
        
        self.check_nodeName()
        
        name = self.get_amplitude_name(f, n=n, m=m)
        
        if n==None and m==None:
            rtn.append("ad {} {} {}".format(name, f, self.nodeName))
        else:
            rtn.append("ad {} {} {} {} {}".format(name, f, n, m, self.nodeName))
            
        return rtn
    
    def _pdtype(self, name, sigtype):
        rtn = []
        
        if sigtype == "pitch":
            rtn.append("pdtype {} y-split".format(name))
        elif sigtype == "yaw":
            rtn.append("pdtype {} x-split".format(name))
        
        return rtn
    
    def add_signal(self, quad=None, sigtype=None):
        """
        Adds a photodiode detector to the kat object at this port. Must
        specify which demodulation quadrature and type of signal to
        detect.
        
        quad: "I" or "Q", Demodoulation quadrature relative to the Output's `phase` value.
        sigtype: "z","pitch" or "yaw", type of signal to detect
        
        Returns: Name of added detector
        
        
        Example:
            REFL_dets = [] # Store the names of each detector
        
            REFL_dets.append( base.IFO.ASC_REFL36A.add_signal('I', "pitch") )
            REFL_dets.append( base.IFO.ASC_REFL36A.add_signal("Q", "pitch") )
            REFL_dets.append( base.IFO.ASC_REFL36B.add_signal("I", "pitch") )
            REFL_dets.append( base.IFO.ASC_REFL36B.add_signal("Q", "pitch") )
        """
        cmds = self.get_signal_cmds(quad=quad, sigtype=sigtype)
        self.__IFO.kat.parse(cmds, addToBlock=self._block)
        return self.get_signal_name(quad, sigtype)
        
    def get_signal_name(self, quad="I", sigtype='z'):
        
        name = self.name
        
        # If we're demodulating add which quadrature we're using
        if self.f is not None: name += "_" + quad
        
        if sigtype == "pitch":
            name += "_P"
        elif sigtype == "yaw":
            name += "_Y"
            
        return name
        
    def get_signal_cmds(self, dof=None, **kwargs):
        """
        Returns the Finesse commands for a detector added to this port's location.
        
        dof: A DOF object which defines the quadrature and signal type for readout
        
        Optionally keyword arguments can be used for manually picking quadrature
        and signal type (overrides any DOF object setting):
        
        quad: "I" or "Q", Demodoulation quadrature relative to the Output's `phase` value.
        sigtype: "z","pitch" or "yaw", type of signal to detect
        
        Returns: List of commands
            cmds = base.IFO.ASC_REFL36B.get_signal_cmds(kat.IFO.CHARD_P)
            # Or
            cmds = base.IFO.ASC_REFL36B.get_signal_cmds(quad="I", sigtype="pitch")
            
            base.parse(cmds)
        """
        
        if dof is not None:
            if dof.quad is not None: quad = dof.quad
            if dof.sigtype is not None: sigtype = dof.sigtype
            
        if "quad" in kwargs:
            quad = kwargs['quad']
        else:
            quad = None
        if "sigtype" in kwargs:
            sigtype = kwargs['sigtype']
            if sigtype is None:
                sigtype = 'z'
        else:
            sigtype = 'z'
            
        if self.f is not None and quad is None:
            raise pkex.BasePyKatException("No quadrature value specified")
        
        rtn = []
        
        self.check_nodeName()
        
        name = self.get_signal_name(quad=quad, sigtype=sigtype)
                
        if self.f==None:
            rtn.append("pd {} {}".format(name, self.nodeName))
        else:
            if quad !="I" and quad != "Q":
                raise pkex.BasePyKatException("quadrature must be 'I' or 'Q'")
                
            phase = self.IQ_phase(quad, self.phase)
            
            rtn.append("pd1 {} {} {} {}".format(name, self.f, phase, self.nodeName))
        
        rtn.extend(self._pdtype(name, sigtype))
        
        return rtn
        
    def IQ_phase(self, quad, phase):
        if phase is None:
            raise pkex.BasePyKatException("Phase cannot be None")
            
        if quad == "Q":
            phase = phase + 90.0
            
            if phase >=360.0 :
                phase -= 360.0
                
        return phase
    
    def add_transfer(self, quad=None, sigtype="z"):
        cmds = self.get_transfer_cmds(quad=quad, sigtype=sigtype)
        self.__IFO.kat.parse(cmds, addToBlock=self._block)
        return self.get_transfer_name(quad, sigtype)
        
    def get_transfer_name(self, quad="I", sigtype="z"):
        name = self.name
        
        # If we're demodulating add which quadrature we're using
        if self.f is not None: name += "_" + quad
        
        if sigtype == "pitch":
            name += "_P"
        elif sigtype == "yaw":
            name += "_Y"
            
        return name + "_TF"
        
    def get_transfer_cmds(self, quad="I", sigtype="z", phase2=None):
        self.check_nodeName()
        
        name = self.get_transfer_name(quad=quad, sigtype=sigtype)
        
        rtn = []
            
        if self.f is None:
            if phase2 is None:
                rtn.append("pd1 {} {} {}".format(name, "$fs", self.nodeName))
            else:
                rtn.append("pd1 {} {} {} {}".format(name, "$fs", phase2, self.nodeName))
        else:
            if quad not in ("I", "Q"):
                raise pkex.BasePyKatException("quadrature must be 'I' or 'Q'")            
                
            phase = self.IQ_phase(quad, self.phase)
            
            if phase2 is None:
                rtn.append("pd2 {} {} {} {} {}".format(name, self.f , phase, "$fs", self.nodeName))
            else:
                rtn.append("pd2 {} {} {} {} {} {}".format(name, self.f, phase, "$fs", phase2, self.nodeName))
        
        rtn.extend(self._pdtype(name, sigtype))
        
        return rtn    
                
    def get_diff_cmds(self, quad='I', sigtype='z'):
        '''
        Generates code for a detector to be used for computing slope of an error signal using the
        diff command in Finesse. Might be unnecessary method, other ways to generate this detector.
        Does the same as port.get_signal_cmds() I think, thus remove when sure. 
        '''
        self.check_nodeName()
        name = self.get_transfer_name(quad=quad, sigtype=sigtype)
        rtn = []
        if self.f is None:
            rtn.append("pd {} {}".format(name, self.nodeName))
        else:
            if quad not in ("I", "Q"):
                raise pkex.BasePyKatException("quadrature must be 'I' or 'Q'")            
                
            phase = self.IQ_phase(quad, self.phase)
            
            rtn.append("pd1 {} {} {} {}".format(name, self.f, phase, self.nodeName))

        rtn.extend(self._pdtype(name, sigtype))
                
        return rtn            
