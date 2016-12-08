from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import math
import copy
import warnings
import cmath
import pykat.components
import pykat.exceptions as pkex
import pykat.external.peakdetect as peak
import matplotlib.pyplot as plt


nsilica=1.44963098985906


class aLIGO(object):
    def __init__(self, _name="default"):
        names = ['default', 'LLO', 'LHO']
        if _name not in names:
            printf("aLIGO name `{}' not recognised, must be 'default', 'LLO' or 'LHO'",_name)

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

def scan_optics(_kat, _optics, _factors, detector, xlimits=[-100, 100], steps=200, minmax='max', debug=False): 
    """
    Scans one or more optics (by changing its tuning) and find the tuning with max/min
    output on a sepcific detectos. Useful for pre-tuning cavities or interferometers.
    Returns the tuning of the maximum (or minimum) and the precision
    Parameters:
    optics: list of names of components to be tuned
    factors: list of scaling factors for the tuning for each optics, first element must be 1.0
    detector: name of detetor to use
    xlimits: limits of scan, defaults to [-100, 100]
    steps: number of steps to use in scan, default is 200
    minmax, string to indicate maximum or minum detection, default 'max'
    Usage:
    scan_optis(kat, "PRM", 1, "pdout")
    scan_optis(kat, ["ETMX", "ETMY", [1, -1], "pdout", minmax='min')
    """
    kat = _kat.deepcopy()
    
    optics=make_list_copy(_optics)
    factors=make_list_copy(_factors)
    if len(optics) != len(factors):
        raise pkex.BasePyKatException("You must provide a factor for each optics")
    if factors[0]!=1.0:
        raise pkex.BasePyKatException("First element in `factors' must be 1.0")

    kat.parseCommands("""
    xaxis {0} phi lin {1} {2} {3}
    """.format(optics[0], xlimits[0], xlimits[1], steps))

    idx=1
    for o in optics[1:]:
        kat.parseCommands("""
        func nx{0} = ({1}) * $x1
        put {2} phi $nx{0}""".format(idx,factors[idx],o))
    
    out = kat.run()
    
    stepsize = out.x[1]-out.x[0]
    
    print("  stepsize (precision) of scan: {0:g}".format(stepsize))

    if debug==True:
        plt.figure()
        plt.plot(out.x,out[detector])

    _max, _min = peak.peakdetect( out[detector],out.x, 1)
    
    if debug==True:
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
