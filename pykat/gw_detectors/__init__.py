import pykat
import pykat.exceptions as pkex
import pykat.external.peakdetect as peak

import numpy as np
import inspect
import math

global nsilica, clight
nsilica = 1.44963098985906
clight = 299792458.0

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
    
    
    
def round_to_n(x, n):
    if not x: return 0
    power = -int(math.floor(math.log10(abs(x)))) + (n - 1)
    factor = (10 ** power)
    return round(x * factor) / factor
    
def vprint(verbose, printstr):
    if verbose:
        print(printstr)
        
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
         
def make_list_copy(_l):
    """
    Utility function, takes a list of strings or single string
    and returns a copy of a list, e.g.
    "string" copy to ["string"]
    ["string1", "string2"] copy to ["string1", "string2"]
    """
    if not isinstance(_l, (list, tuple)):
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
                if _ in self.__DOFs:
                    dofs.append(self.__DOFs[_])
                else:
                    raise pkex.BasePyKatException("Could not find DOF called `%s`. Possible DOF options: %s" % (_, str(list(self.__DOFs.keys()))))
            else:
                raise pkex.BasePyKatException("'%s' not possible DOF options: %s" % (_, str(list(self.__DOFs.keys()))))
        
        return dofs
    
    def scan_DOF_cmds(self, DOF, xlimits=[-100, 100], steps=200, relative=False):
        return scan_optics_string(DOF.optics, 
                                  DOF.factors,
                                  "scan",
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
                                  
    def scan_DOF(self, DOF, xlimits=[-100, 100], steps=200, relative=False): 
        kat = self.kat.deepcopy()
        kat.parseCommands(self.scan_DOF_cmds(DOF, xlimits=xlimits, steps=steps, relative=relative))
        kat.parseCommands(DOF.signal())
        
        return kat.run()

    def scan_optics(self, _optics, _factors, xlimits=[-100, 100], steps=200,relative=False): 
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
        kat = self.kat.deepcopy()
    
        optics=make_list_copy(_optics)
        factors=make_list_copy(_factors)
    
        kat.parseCommands(self.scan_optic_cmds(_optics, _factors, xlimits=xlimits, steps=steps, relative=relative))
    
        return kat.run(cmd_args="-cr=on")

    def optical_gain(self, DOF_sig, DOF_det, f=10.0):
        kat = self.kat.deepcopy()
        kat.removeBlock('locks', False)
         
        _fsigStr = DOF_sig.fsig("sig1", fsig=f)
        _detStr = DOF_det.transfer(fsig=f)
        _detName = DOF_det.transfer_name()
        
        kat.parseCommands(_fsigStr)
        kat.parseCommands(_detStr)
        kat.noxaxis = True
        kat.parseCommands("yaxis lin abs:deg")
        
        out = kat.run()
        
        return float(np.real(out[_detName]))
        
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

    def apply_tuning(self, phi, add=False):
        for idx, o in enumerate(self.optics):
            if add:
                self.__IFO.kat.components[o].phi += phi * self.factors[idx]
            else:
                self.__IFO.kat.components[o].phi = phi * self.factors[idx]
            
    def signal(self):
        return self.port.signal(self.quad, sigtype=self.sigtype)
        
    def signal_name(self):
        return self.port.signal_name(self.quad, sigtype=self.sigtype)

    def transfer(self, fsig, phase2=None):
        return self.port.transfer(self.quad, fsig=fsig, phase2=phase2, sigtype=self.sigtype)
        
    def transfer_name(self):
        return self.port.transfer_name(self.quad)

    def fsig(self, _fsigName, fsig=1.0):
        _fsigStr= ""
        
        for idx, o in enumerate(self.optics):
            phase = 0.0
            if self.factors[idx] == -1:
                phase = 180.0
            _fsigStr = "\n".join([_fsigStr, "fsig {} {} {} {} ".format(_fsigName, o, fsig, phase)])
            
        return _fsigStr

class Port(object):
    """
    Defining an output port for the interferometer, can be either a
    pd or a pd1 detector (for error signal generation).
    """
    def __init__(self, IFO, _portName, _nodeNames, f=None, phase=None):
        self.__IFO = IFO
        self.portName = _portName
        self.nodeNames = make_list_copy(_nodeNames)
        self.f=f            # demodulation frequency, float
        self.phase = phase  # demodulation frequency for I quadrature, float
        self.name = self.portName            
    
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

    def amplitude_name(self, f, n=None, m=None, sigtype="z"):
        name = self.name + "_ad"
        return name
    
    def amplitude(self, f, n=None, m=None, sigtype="z"):
        self.check_nodeName()
        
        name = self.amplitude_name(f, n=n, m=m, sigtype=sigtype)
        
        if n==None and m==None:
            return "ad {} {} {}".format(name, f, self.nodeName)
        else:
            return "ad {} {} {} {} {}".format(name, f, n, m, self.nodeName)
    
    def signal_name(self, quad="I", sigtype="z"):
        name = self.name
        
        if self.f!=None:
            name = self.name+"_"+quad
            
        return name
    
    def signal(self, quad="I", sigtype="z"):
        self.check_nodeName()
        
        name = self.signal_name(quad=quad, sigtype=sigtype)
        
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
        
    def transfer_name(self, quad="I"):
        name = self.name
        
        if self.f!=None:
            name = self.name+"_"+quad
            
        return name
        
    def transfer(self, quad="I", fsig=1.0, phase2=None, sigtype="z"):
        self.check_nodeName()
        
        name = self.transfer_name(quad=quad)
        
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
                
                