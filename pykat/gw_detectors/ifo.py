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

        set_lengths(self)
            
        # useful ports
        self.POP_f1  = port("POP",   "nPOP",  "f1", phase=101)
        self.POP_f2  = port("POP",   "nPOP",  "f2", phase=13)
        self.REFL_f1 = port("REFL",  "nREFL", "f1", phase=101)
        self.REFL_f2 = port("REFL",  "nREFL", "f2", phase=14)
        self.AS_DC   = port("AS_DC", "nSRM2")
        self.POW_BS  = port("PowBS", "nPRBS*")
        self.POW_X   = port("PowX",  "nITMX2")
        self.POW_Y   = port("PowY",  "nITMY2")

        # pretune DOF
        self.preARMX =  DOF("ARMX", self.POW_X,   "", "ETMX", 1, 1.0)
        self.preARMY =  DOF("ARMY", self.POW_Y,   "", "ETMY", 1, 1.0)
        self.preMICH =  DOF("AS"  , self.AS_DC,   "", "BS",   1, 6.0)
        self.prePRCL =  DOF("PRCL", self.POW_BS,  "", "PRM",  1, 10.0)
        
        # control scheme as in C. Bond Table C.1,  TODO complete citation
        self.PRCL =  DOF("PRCL", self.POP_f1,  "I", "PRM", 1, 100.0)
        self.MICH =  DOF("MICH", self.POP_f2,  "Q", "BS",  1, 100.0)
        self.CARM =  DOF("CARM", self.REFL_f1, "I", ["ETMX", "ETMY"], [1, 1], 1.5)
        self.DARM =  DOF("DARM", self.AS_DC,   "",  ["ETMX", "ETMY"], [1,-1], 1.0)
        self.SRCL =  DOF("SRCL", self.REFL_f2, "I", "SRM", 1, 1e2)

        # Defining a dicotionary for the main mirror positions (tunings)
        self.tunings = {}
        self.tunings = self.get_tunings()

    def set_lengths(self, kat, verbose=False):
        # distances between HR surfaces:
        self.lpr = kat.lp1.L + kat.lp2.L + kat.lp3.L
        self.lx = kat.lx.L + kat.BSsub1.L * kat.BSsub1.n + kat.ITMXsub.L * kat.ITMXsub.n
        self.ly = kat.ly.L + kat.ITMYsub.L * kat.ITMYsub.n
        self.lsr = kat.ls1.L + kat.ls2.L + kat.ls3.L + kat.BSsub2.L * kat.BSsub2.n
        # resulting combined distances (single, not roundtrip)
        self.lMI =  0.5 * (self.lx + self.ly)
        self.lPRC = self.lpr + self.lMI
        self.lSRC = self.lsr + self.lMI
        self.lSchnupp = self.lx - self.ly
        if verbose:
            print("-- small MI lengths")
            print(" lx = {}m, ly = {}m".format(self.lx, self.ly))
            print(" lpr = {}m, lsr = {}m".format(self.lpr, self.lsr))
            print(" lMI = {}m, lSchnupp = {}m".format(self.lMI, self.lSchnupp))
            print(" lpr = {}m".format(self.lpr))
            print(" lpr = {}m".format(self.lpr))
        
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

    def pretune(self, _kat, pretune_precision=1.0e-4):
        print("-- pretuning interferometer to precision {0:2g} deg = {1:2g} m".format(pretune_precision, pretune_precision*_kat.lambda0/360.0))
        kat=_kat.deepcopy()
        _kat.ITMX.phi=0.0
        _kat.ITMY.phi=0.0
        #print("  removing components, commands and blocks")
        #remove_components(kat, ["mod1", "lmod2", "mod2", "lmod3"], component_in="lmod1");
    
        print("  scanning X arm")
        kat1 = kat.deepcopy()
        make_transparent(kat1,["PRM","SRM"])
        make_transparent(kat1,["ITMY", "ETMY"])
        kat1.BS.setRTL(0.0,1.0,0.0) # set BS refl. for X arm
        phi, precision = self.scan_to_precision(kat1, self.preARMX, pretune_precision)
        phi=round(phi*pretune_precision)/pretune_precision
        print("  found max/min at: {} (precision = {:2g})".format(phi, precision))
        _kat.ETMX.phi=kat.ETMX.phi=phi
    
        print("  scanning Y arm")
        kat1 = kat.deepcopy()
        make_transparent(kat1,["PRM","SRM"])
        make_transparent(kat1,["ITMX", "ETMX"])
        kat1.BS.setRTL(1.0,0.0,0.0) # set BS refl. for Y arm
        phi, precision = self.scan_to_precision(kat1, self.preARMY, pretune_precision)
        phi=round(phi*pretune_precision)/pretune_precision
        print("  found max/min at: {} (precision = {:2g})".format(phi, precision))
        _kat.ETMX.phi=kat.ETMX.phi=phi
    
        print("  scanning MICH")
        kat1=kat.deepcopy()
        make_transparent(kat1,["PRM","SRM"])
        phi, precision = self.scan_to_precision(kat1, self.preMICH, pretune_precision, minmax="min", precision=30.0)
        phi=round(phi*pretune_precision)/pretune_precision
        print("  found max/min at: {} (precision = {:2g})".format(phi, precision))
        _kat.BS.phi=kat.BS.phi=phi

        print("  scanning PRCL")
        kat1=kat.deepcopy()
        make_transparent(kat1,["SRM"])
        phi, precision = self.scan_to_precision(kat1, self.prePRCL, pretune_precision, precision=30.0)
        phi=round(phi*pretune_precision)/pretune_precision
        print("  found max/min at: {} (precision = {:2g})".format(phi, precision))
        _kat.PRM.phi=kat.PRM.phi=phi

        print("  scanning SRCL")
        kat1=kat.deepcopy()
        """
        phi = -90 # start near expected RSE point
        kat1.ETMY.phi += 1e-6 # applying some random DC offset
        kat1.ETMX.phi -= 1e-6
        TFstr, name = self.DARM.port.transfer(kat1)
        kat1.parseCommands(TFstr)
        kat1.parseCommands(self.DARM.fsig("sig1"))
        kat1.parseCommands("xaxis SRM phi lin 0 10 200")
        precision = 60.0
        while precision>pretune_precision:
            kat1.xaxis.limits=[phi-1.5*precision, phi+1.5*precision]
            out = kat1.run()
            phi, precision = find_peak(out, self.DARM.port.name, minmax="max")
        """
        phi, precision = self.scan_to_precision(kat1, self.prePRCL, pretune_precision, precision=30.0, phi=0)
        phi=round(phi*pretune_precision)/pretune_precision - 90.0
        print("  found max/min at: {} (precision = {:2g})".format(phi, precision))
        _kat.SRM.phi=kat.SRM.phi=phi
        
                
        tuning = self.get_tunings(_kat)
        return tuning
        # done
        
    def scan_to_precision(self, kat, DOF, pretune_precision, minmax="max", phi=0.0, precision=60.0):
        while precision>pretune_precision*DOF.scale:
            out = scan_DOF(kat, DOF, xlimits = [phi-1.5*precision, phi+1.5*precision])
            phi, precision = find_peak(out, DOF.port.portName, minmax=minmax)
            #print("** phi= {}".format(phi))
        return phi, precision


    def power_ratios(self, _kat):
        kat = _kat.deepcopy()
        kat.verbose = False
        kat.noxaxis = True
    
        ports = [self.POW_X, self.POW_Y, self.AS_DC, self.POW_BS]
        _detStr = ""
        for p in ports:
            _sigStr, name = p.signal(kat)
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
            ax.semilogy(out.x,out[d.port.name])
            ax.set_xlim([np.min(out.x), np.max(out.x)])
            ax.set_xlabel("phi [deg] {}".format(d.optics[0]))
            ax.set_ylabel('{} [W] '.format(d.port.name))
            ax.grid()
        plt.tight_layout()
        plt.show(block=0)

    def plot_error_signals(self, _kat, xlimits=[-1,1]):
        kat = _kat.deepcopy()
        kat.verbose = False
        kat.noxaxis = True
        dofs = [self.CARM, self.DARM, self.PRCL, self.MICH]
        idx=1
        fig = plt.figure()
        for d in dofs:
            ax = fig.add_subplot(2,2,idx)
            idx+=1
            out = scan_DOF(kat, d, xlimits = np.multiply(d.scale,xlimits), relative = True)
            ax.plot(out.x,out[d.port.name])
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
        _sigStr, name = self.AS_DC.signal(kat)
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
        return round(out[0],6)
            
class port(object):
    """
    Defining an output port for the interferometer, can be either a
    pd or a pd1 detector (for error signal generation).
    """
    def __init__(self, _portName, _nodeNames, f=None, phase=None):
        self.portName = _portName
        self.nodeNames = make_list_copy(_nodeNames)
        self.f=f            # demodulation frequency, string "f1", "f2" or "f3"
        self.phase = phase  # demodulation frequency for I quadrature, float
        if f!=None:
            self.name = self.portName+"_"+str(self.f)
        else:
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

    def signal(self, kat, quad="I", sigtype="z"):
        self.check_nodeName(kat)
        name = self.name
        if sigtype != "z":
                raise pkex.BasePyKatException("alignment signals are not implemented yet")            
        if self.f==None:
            return "pd {} {}".format(self.name, self.nodeName), name
        else:
            if quad !="I" and quad != "Q":
                raise pkex.BasePyKatException("quadrature must be 'I' or 'Q'")            
            name = self.name+"_"+quad
            phase = self.IQ_phase(quad, self.phase)
            return "pd1 {} ${} {} {}".format(self.name, "f1", self.phase, self.nodeName), name
        
    def IQ_phase(self, quad, phase):
        if quad== "Q":
            phase = phase + 90.0
            if phase >=360.0 :
                phase -= 360.0
        return phase


    def amplitude(self, kat, f, n=None, m=None, sigtype="z"):
        self.check_nodeName(kat)
        name = self.name + "_ad"
        if n==None and m==None:
            return "ad {} {} {}".format(self.name, f, self.nodeName), name
        else:
            return "ad {} {} {} {} {}".format(self.name, f, n, m, self.nodeName), name
        
    def transfer(self, kat, quad="I", fsig=1.0, phase2=None, sigtype="z"):
        name = self.name
        self.check_nodeName(kat)
        if sigtype!="z":
                raise pkex.BasePyKatException("alignment signals are not implemented yet")            
        if self.f==None:
            if phase2 == None:
                return "pd1 {} {} {}".format(name, fsig, self.nodeName), name
            else:
                return "pd1 {} {} {} {}".format(name, fsig, phase2, self.nodeName), name

        else:
            if quad !="I" and quad != "Q":
                raise pkex.BasePyKatException("quadrature must be 'I' or 'Q'")            
            name=self.name+"_"+quad
            phase = IQ_phase(quad, self.phase)
            if phase2 == None:
                return "pd2 {} ${} {} {} {}".format(self.name, "f1", self.phase, fsig, self.nodeName), name
            else:
                return "pd2 {} ${} {} {} {} {}".format(self.name, "f1", self.phase, fsig, phase2, self.nodeName), name
                                
class DOF(object):
    """
    Defining a degree of freedom for the interferometer, includes the
    objects and how to move them, and the default output port to read
    out the DOF signal.
    """
    def __init__(self, _DOFName, _port, _quad, _optics, _factors, _scale):
        self.name = _DOFName
        self.port = _port
        self.quad = _quad
        self.optics=make_list_copy(_optics)
        self.factors=make_list_copy(_factors)
        # scaling factor, to compensate for lower sensitivity compared
        # to DARM (in tuning plots for example)
        # Thus DARM has a scale of 1, all other DOFs a scale >1
        self.scale = _scale
        
    def fsig(self, _fsigName):
        _fsigStr= ""
        for idx, o in enumerate(self.optics):
            phase = 0.0
            if self.factors[idx] == -1:
                phase = 180.0
            _fsigStr = "\n".join([_fsigStr, "fsig {} {} 1 {} ".format(_fsigName, o, phase)])
        return _fsigStr

        
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
    sigStr, name = DOF.port.signal(kat)
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
