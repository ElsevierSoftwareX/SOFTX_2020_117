import pykat
import pykat.exceptions as pkex
from pykat import isContainer
import numpy as np
import inspect
import math
import six
from copy import deepcopy

from pandas import DataFrame

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

            ax.plot((theta,theta), (r_lim[0], r), label=re.sub("[\(\[].*?[\)\]]", "", _).strip(), lw=2)

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
    
def vprint(verbose, printstr):
    if verbose:
        print(printstr)
        
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
        plt.plot(out.x,out[detector])
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


def scan_optics_string(_optics, _factors, _varName, target="phi", linlog="lin", xlimits=[-100, 100], steps=200, axis=1,relative=False):
    optics=make_list_copy(_optics)
    factors=make_list_copy(_factors)
    
    if len(optics) != len(factors):
        raise pkex.BasePyKatException("you must provide a factor for each optics")

    if linlog not in ["lin", "log"]: 
        raise pkex.BasePyKatException("linlog must be 'lin' or 'log'")
        
    _tuneStr  = "var {} 0\n".format(_varName)
    
    if axis==1:
        _tuneStr += "xaxis {} {} {} {} {} {}".format(_varName, target, linlog, xlimits[0], xlimits[1], steps)
    elif (axis==2 or axis==3): 
        _tuneStr += "x{}axis {} {} {} {} {} {}".format(axis, _varName, target, linlog, xlimits[0], xlimits[1], steps)
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
                
        _putStr = "\n".join([_putStr, "{} {} {} {}{}".format(_putCmd, o, target, _xStr, axis)])            
        
    _tuneStr += _putStr
        
    return _tuneStr
    
    
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
                              
def scan_DOF(kat, DOF, xlimits=[-100, 100], steps=200, relative=False): 
    kat = kat.deepcopy()
    kat.parse(scan_DOF_cmds(DOF, xlimits=xlimits, steps=steps, relative=relative))
    
    if DOF.port is not None:
        kat.parse(DOF.signal())
    
    return kat.run(cmd_args=["-cr=on"])

def scan_optics(kat, _optics, _factors, target="phi", xlimits=[-100, 100], steps=200,relative=False): 
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

    return kat.run(cmd_args="-cr=on")

def optical_gain(DOF_sig, DOF_det, f=1.0):
    """
    Returns W/rad for length sensing. Will throw an exception if it can't convert the
    units into W/rad.
    """
    
    kat = DOF_sig.kat.deepcopy()
    kat.removeBlock('locks', False)
     
    _fsigStr = DOF_sig.fsig("sig1", fsig=f)
    _detStr  = DOF_det.transfer()
    _detName = DOF_det.transfer_name()
    
    kat.parse(_fsigStr)
    kat.parse(_detStr)
    kat.noxaxis = True
    kat.parse("yaxis lin abs:deg")
    
    out = kat.run()
    
    if DOF_sig.sigtype == "phase":
        return float(np.real(out[_detName])) # W/rad
    elif DOF_sig.sigtype == "z":
        k = 2*np.pi/kat.lambda0
        return float(np.real(out[_detName])) / k # W/(m*k) = W/rad
    else:
        raise pkex.BasePyKatException("Not handling requested sigtype for unit conversion")

def diff_DOF(DOF, target, deriv_h=1e-12):
    """
    Returns commands to differentiate with respect to the DOF motion.
    
    This is typically used to find the slope of error signals at DC.
    """
    
    if target == "z":
        # As we can't target the z position of the component directly
        # we need to target phi and scale things
        pass
        
    rtn = ("var x 0\n"
            "diff x re\n"
            "deriv_h {deriv_h}\n" 
            "set _dx x re\n"
            "func DX = ({scaling}) * $_dx\n"
            "noplot DX\n"
            "func mDX = (-1 * {scaling}) * $_dx\n"
            "noplot mDX\n").format(deriv_h=deriv_h, scaling=scaling)
    
    for o,f in zip(DOF.optics, DOF.factors):
        if f == 1:
            rtn += "put %s %s $DX\n" % (o,target)
        elif f == -1:
            rtn += "put %s %s $mDX\n" % (o,target)
        else:
            raise pkex.BasePyKatException("Factor can only be -1 or 1 currently")
    
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
    
def optimise_demod_phase(_kat, DOF, ports, debug=False):
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
            _ports.append(kat.IFO.Ports[port])
        else:
            _ports.append(kat.IFO.Ports[port.name])
            
        if _ports[-1].f is None:
            raise pkex.BasePyKatException("port %s cannot have its demodulation phase optimised as it isn't demodulated" % port.name)
    
    ports = _ports
    
    pd_detectors = []
        
    # Add in the signals
    for port in ports:
        pd_detectors.append(kat.IFO.Ports[port.name].add_signal("I", DOF.sigtype))
    
    if debug:
        print("Optimising pds: %s" % pd_detectors)
        
    kat.removeBlock("locks", False)
    kat.removeBlock("powers", False)
    
    kat.parse( aligo.diff_DOF(DOF, DOF.sigtype), addToBlock="OPTIMISE")
    kat.parse( aligo.scan_demod_phase_cmds(pd_detectors), addToBlock="OPTIMISE")
    
    # Analyitcally we can find the phase which gives a maxmium
    # by solving for A and B in:
    #   f(x) = Acos(x+B)
    # assuming we take the two data points x = {0, pi/2} we find
    #   B = arctan(y2/y1)
    kat.xaxis.limits = (0, 90)
    kat.xaxis.steps = 1
    
    if debug:
        print(kat & "OPTIMISE")
        
    out = kat.run()
    
    if debug:
        print(out.y)
        
    rtn = []
    
    for pd, port in zip(pd_detectors, ports):
        
        y1, y2 = out[pd]
        x = np.deg2rad(out.x)
        R = np.sqrt(y1**2 + y2**2)
        phi = np.rad2deg(np.arctan2(y2,y1))
        
        # All in I quadrature so no
        _kat.IFO.Ports[port.name].phase = phi
        
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
    
    cavs = []
    qxs = []
    qys = []

    for cav in _kat.getAll(pykat.commands.cavity):
        # Switch off all cavities
        for _ in _kat.getAll(pykat.commands.cavity):
            _.enabled = False
        # Then select one at a time
        cav.enabled = True

        out = _kat.run()

        cavs.append(cav.name)
        qxs.append(pykat.BeamParam(q=out['qx']))
        qys.append(pykat.BeamParam(q=out['qy']))

    mmx = DataFrame(index=cavs, columns=cavs)
    mmy = DataFrame(index=cavs, columns=cavs)
    
    for c1, q1x, q1y in zip(cavs, qxs, qys):
        for c2, q2x, q2y in zip(cavs, qxs, qys):
            mmx[c1][c2] = pykat.BeamParam.overlap(q1x, q2x)
            mmy[c1][c2] = pykat.BeamParam.overlap(q1y, q2y)
            
    return 1-mmx, 1-mmy, list(zip(cavs, qxs, qys))

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
        
        See `scan_DOF` for keyword arguments options.
        """
        return scan_DOF(self.__IFO.kat, self, **kwargs)
        
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

class Port(object):
    """
    This object defines a location in an interferometer where detectors are places and demodulated at a particular
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
    def __init__(self, IFO, _portName, _nodeNames, f=None, phase=0, block=None):
        self.__IFO = IFO
        self.portName = _portName
        self.nodeNames = make_list_copy(_nodeNames)
        self.f = f            # demodulation frequency, float
        self.phase = phase    # demodulation phase for I quadrature, float
        self.name = self.portName            
        self._block = block

    @property
    def kat(self):
        """For referencing the kat object this DOF is associated with"""
        return self.__IFO.kat
          
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

    def get_amplitude_name(self):
        return self.name + "_ad"
    
    def get_amplitude_cmds(self, f, n=None, m=None):
        rtn = []
        
        self.check_nodeName()
        
        name = self.amplitude_name(f, n=n, m=m)
        
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
        
        quad: "I" or "Q", Demodoulation quadrature relative to the Port's `phase` value.
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
        
        quad: "I" or "Q", Demodoulation quadrature relative to the Port's `phase` value.
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
            
        if "quad"    in kwargs: quad    = kwargs['quad']
        if "sigtype" in kwargs: sigtype = kwargs['sigtype']
            
        if self.f is not None and quad is None: raise pkex.BasePyKatException("No quadrature value specified")
        if sigtype is None: sigtype = "z"
        
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
            
        return name
        
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
                