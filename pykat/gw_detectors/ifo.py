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

nsilica=1.44963098985906

class aLIGO(object):
    """
    Object storing the Finesse input file and some auxilliary information
    about Avanced LIGO interferometers.
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

        self.POP_f1  = port("POP",  "nPOP",  "f1", phase=101)
        self.POP_f2  = port("POP",  "nPOP",  "f2", phase=13)
        self.REFL_f1 = port("REFL", "nREFL", "f1", phase=101)
        self.REFL_f2 = port("REFL", "nREFL", "f2", phase=14)
        self.AS_DC   = port("AS",   "nSRM")

        # control scheme as in C. Bonf Table C.1 TODO complete citation
        self.PRCL =  DOF("PRCL", self.POP_f1,  "I", "PRM", 1)
        self.MICH =  DOF("MICH", self.POP_f2,  "Q", "BS",  1)
        self.CARM =  DOF("CARM", self.REFL_f1, "I", ["ETMX", "ETMY"], [1, 1])
        self.DARM =  DOF("DARM", self.AS_DC,   "",  ["ETMX", "ETMY"], [1,-1])
        self.SRCL =  DOF("SRCL", self.REFL_f2, "I", "SRM", 1)
        
        

class port(object):
    """
    Defining an output port for the interferometer, can be either a
    pd or a pd1 detector (for error signal generation).
    """
    def __init__(self, _portName, _nodeName, f=None, phase=None):
        self.portName = _portName
        self.nodeName = _nodeName
        self.f=f            # demodulation frequency, string "f1", "f2" or "f3"
        self.phase = phase  # demodulation frequency for I quadrature, float
        if f!=None:
            self.name = self.portName+"_"+str(self.f)
        
    def errsig(quad="I", sigtype="z"):
        if sigtype != "z":
                raise pkex.BasePyKatException("alignment signals are not implemented yet")            
        if self.f==None:
            return "pd {} {}".format(name, self.nodeName)
        else:
            if quad !="I" and quad != "Q":
                raise pkex.BasePyKatException("quadrature must be 'I' or 'Q'")            
            name = self.name+"_"+quad
            phase = IQ_phase(quad, self.phase)
            return "pd1 {} ${} {} {}".format(name, f1, phase, self.nodeName), name
        
    def IQ_phase(quad, phase):
        if quad== "Q":
            phase = phase + 90.0
            if phase >=360.0 :
                phase -= 360.0
        return phase

    def amplitude(f, n=None, m=None, sigtype="z"):
        name = self.name + "_ad"
        if n==None and m==None:
            return "ad {} {} {}".format(name, f, self.nodeName)
        else:
            return "ad {} {} {} {} {}".format(name, f, n, m, self.nodeName), name
        
    def transfer(quad="I", fsig=1, phase2=None, sigtype="z"):
        if sigtype!="z":
                raise pkex.BasePyKatException("alignment signals are not implemented yet")            
        if self.f==None:
            return "pd1 {} {} {} {}".format(name, fsig, phase2, self.nodeName), name
        else:
            if quad !="I" and quad != "Q":
                raise pkex.BasePyKatException("quadrature must be 'I' or 'Q'")            
            name=self.name+"_"+quad
            phase = IQ_phase(quad, self.phase)
            if phase2 == None:
                return "pd2 {} ${} {} {} {}".format(name, f1, phase, fsig, self.nodeName), name
            else:
                return "pd2 {} ${} {} {} {} {}".format(name, f1, phase, fsig, phase2, self.nodeName), name
                    
        
class DOF(object):
    """
    Defining a degree of freedom for the interferometer, includes the
    objects and how to move them, and the default output port to read
    out the DOF signal.
    """
    def __init__(self, _DOFName, _port, _quad, _optics, _factors):
        self.DOFName = _DOFName
        self.port = _port
        self.quad = _quad
        self.optics=make_list_copy(_optics)
        self.factors=make_list_copy(_factors)
        
        def tune(_varName, linlog="lin", xlimits=[-100, 100], steps=200, axis=1):
            if linlog not in ["lin", "log"]: 
                raise pkex.BasePyKatException("linlog must be 'lin' or 'log'")
            _tuneStr  = "var {} 0\n".format(_varName)
            if axis==1:
                _tuneStr += "xaxis {} phi {} {} {} {}\n".format(_varName, linlog, xlimits[0], xlimits[1], steps)
            elif (axis==2 or axis==3): 
                _tuneStr += "x{}axis {} phi {} {} {} {}\n".format(axis, _varName, linlog, xlimits[0], xlimits[1], steps)
            else:
                raise pkex.BasePyKatException("axis must be 1, 2 or 3")
            _putStr = ""
            for idx, o in enumerate(self.optics):
                if self.factors[idx] == 1:
                    _xStr="$x"
                elif self.factors[idx] == -1:
                    _xStr="$mx"
                else:
                    raise pkex.BasePyKatException("optics factors must be 1 or -1")
                _putStr = "\n".join([_putStr, "put* {} phi {}{}".format(o, _xStr, axis)])
            _tuneStr += _putStr
            return _tuneStr
                
        def fsig(_fsigName):
            _fsigStr= ""
            for idx, o in enumerate(self.optics):
                _fsigStr = "\n".join([_fsigStr, "fsig {} {} {}".format(_fsigName, o, self.factors[idx])])

            return _fsigStr

    
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
            print('  made {} transparent'.format(o))
    if len(components) != 0:
        raise pkex.BasePyKatException("Cannot find component {}".format(components))
    return kat

def magic_length(): # TODO: improve and reference this!
    # Cavity lengths
    #const Lprc 57.656
    #const Lsrc 56.008
    #const Lschnupp 0.08

    # Individual lengths
    # PRC
    #const Lpr1 16.6107
    #const Lpr2 16.1647
    #const Lpr3 19.5381 

    # SRC
    #const Lsr1 15.7586
    #const Lsr2 15.4435
    #const Lsr3 19.3661

    Laver = 57.656 - 16.6107 - 16.1647 - 19.5381
    # x length between BS and ITM
    Lmx = Laver + 0.5*0.08 - 0.06873 * 1.44963098985906 - 0.2*1.44963098985906
    Lmy = Laver - 0.2*1.44963098985906 - 0.5*0.08
    #func Lasrc = $Laver + $BSthickness * $nsilica
    #func Lsr3 = $Lsrc - $Lsr1 - $Lsr2 - $Lasrc
    #put ls3 L $Lsr3

    return Lmx, Lmy

def reconnect_nodes(kat, component1, idx1, node_name):
    c_string = component1.getFinesseText()
    c = c_string[0].split()
    new_string = " ".join(c[:-2])
    nodes = {}
    nodes[0] = c[-2]
    nodes[1] = c[-1]
    nodes[idx1]=node_name
    new_string = new_string + " " + nodes[0] + " " + nodes[1]
    print(" new string ='{}'".format(new_string))
    kat.parseCommands(new_string)

def remove_commands(kat, _commands):
    commands=make_list_copy(_commands)
    # removing commands
    for o in kat.commands.values():
        if o.name in commands:
            o.remove()
            commands = [c for c in commands if c != o.name]
            print('  {} removed'.format(o))
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
            print('  {} removed'.format(o))
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
    if len(optics) != len(factors):
        raise pkex.BasePyKatException("You must provide a factor for each optics")
    if factors[0]!=1.0:
        raise pkex.BasePyKatException("First element in `factors' must be 1.0")

    if relative == False:
        kat.parseCommands("""
        xaxis {0} phi lin {1} {2} {3}
        """.format(optics[0], xlimits[0], xlimits[1], steps))
    elif relative == True:
        kat.parseCommands("""
        xaxis* {0} phi lin {1} {2} {3}
        """.format(optics[0], xlimits[0], xlimits[1], steps))

    idx=1
    for o in optics[1:]:
        kat.parseCommands("""
        func nx{0} = ({1}) * $x1
        put {2} phi $nx{0}""".format(idx,factors[idx],o))
    
    out = kat.run()
    return out

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
        plt.xlabel('{0} tuning [deg]'.format(optics[0]))
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
