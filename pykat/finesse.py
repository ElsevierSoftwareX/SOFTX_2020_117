# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 09:56:53 2013

PyKat - Python interface and wrapper for FINESSE
Copyright (C) 2013 Daniel David Brown

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

Contact at ddb@star.sr.bham.ac.uk

@author: Daniel Brown
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import subprocess
import tempfile
import numpy as np
import datetime
import time
import pickle
import pykat
import warnings
import re
import math       
import itertools
import ctypes
import ctypes.util
import collections
import re
import copy


try:
    # Python 2
    from itertools import izip_longest
except ImportError:
    # Python 3
    from itertools import zip_longest as izip_longest
"""
try:
    from future_builtins import zip_longest
except ImportError: # not 2.6+ or is 3.x
    try:
        from itertools import izip_longest as zip_longest # < 2.5 or 3.x
    except ImportError:
        print("boom")
        pass
"""

from math import erfc, pi
from collections import namedtuple, OrderedDict

from pykat.node_network import NodeNetwork
from pykat.detectors import BaseDetector as Detector
from pykat.components import Component
from pykat.commands import Command, xaxis
from pykat.SIfloat import *
from pykat.param import Param, AttrParam

import pykat.external.six as six

import pykat.exceptions as pkex

from pykat import USE_GUI, HAS_OPTIVIS, NoGUIException

    
if HAS_OPTIVIS:
    from optivis.bench.labels import Label as optivis_label
    from optivis.geometry import Coordinates as optivis_coord
    import PyQt4

if USE_GUI:
    from pykat.gui.gui import pyKatGUI
    from PyQt4.QtCore import QCoreApplication
    from PyQt4.QtGui import QApplication

from multiprocessing import Process, Manager

NO_BLOCK = "NO_BLOCK"
pykat_web = "www.gwoptics.org/pykat"

# containers used in the trace routine
space_trace = namedtuple("space_trace", ['gouyx','gouyy'])
node_trace = namedtuple("node_trace", ['qx','qy'])
cav_trace = namedtuple("cav_trace", ['isStable','gx','gy','qx','qy','finesse','loss','length','FSR','FWHM','pole'])
         
lkat_location = ctypes.util.find_library("kat")
                                     
def f__lkat_process(callback, cmd, kwargs):
    """
    """
    
    if lkat_location == None:
        raise RuntimeError("Could not find shared library 'libkat', please install to a system location or copy to the same directory as this script")
        
    lkat = ctypes.PyDLL(lkat_location)

    try:
        lkat._pykat_preInit() # must always be called, sets up
                        # exception handling and such no simulation
                        # specifc code here

        # reads in the kat.ini and setups up other parts
        lkat._pykat_init()
        lkat._pykat_setup(cmd)
    
        callback(lkat, **kwargs)
    
    except Exception as ex: 
        print ("Exception caught in python: ", ex.message)
    finally:
        # This should always be called no matter what
        lkat._pykat_finish(0)


def f__lkat_trace_callback(lkat, trace_info, getCavities, getNodes, getSpaces):
    """
    lkat callback for computing the beam traces through a setup.
    Returns a dictionary of nodes, spaces and cavities and the
    various outputs of the tracing algorithm.
    """
    import pylibkat

    # first we need to get a handle on the internals of Finesse
    inter = pylibkat.interferometer.in_dll(lkat, "inter")

    lkat._pykat_step()

    if getNodes:
        for n in range(0, inter.num_nodes):
            node = inter.node_list[n]

            node_info = node_trace(
                                    qx = complex(node.qx.re, node.qx.im),
                                    qy = complex(node.qy.re, node.qy.im)
                                    )

            trace_info[node.name] = node_info

    if getCavities:
        for c in range(0, inter.num_cavities):
            cav = inter.cavity_list[c]

            cav_info = cav_trace(
                                isStable = (cav.stable == 1),
                                gx = cav.stability_x,
                                gy = cav.stability_y,
                                qx = complex(cav.qx.re, cav.qx.im),
                                qy = complex(cav.qy.re, cav.qy.im),
                                finesse = cav.finesse,
                                FSR = cav.FSR,
                                FWHM = cav.FWHM,
                                loss = cav.loss,
                                length = cav.length,
                                pole = cav.pole
                                )

            trace_info[cav.name] = cav_info

    if getSpaces:
        for s in range(0, inter.num_spaces):
            space = inter.space_list[s]

            trace_info[space.name] = space_trace(gouyx = space.gouy_x,
                                                 gouyy = space.gouy_y)
                     
                                             
def GUILength(L):
    """
    Should scale the lengths in some way to handle km and mm for time being
    """
    return L # * ( 40 * erfc(L/400.0) + 0.01)
            
class katRun(object):
    def __init__(self):
        self.runtime = None
        self.StartDateTime = datetime.datetime.now()
        self.x = None
        self.y = None
        self.xlabel = None
        self.ylabels = None
        self.katScript = None
        self.katVersion = None
        self.yaxis = None
        
    def plot(self, logy=False):
        import pylab
        
        pylab.plot(self.x, self.y)
        pylab.legend(self.ylabels,0)
        pylab.xlabel(self.xlabel)
        pylab.show()
        
    def savekatRun(self, filename):
        with open(filename,'w') as outfile:
            pickle.dump(self, outfile)
    
    @staticmethod
    def loadKatRun(filename):
        with open(filename,'r') as infile:
            return pickle.load(infile)
    
    def get(self, value): return self[value]
    
    def __getitem__(self, value):
        idx = [i for i in range(len(self.ylabels)) if self.ylabels[i].split()[0] == str(value)]
        out = None
        
        if len(idx) > 0:
            #out = self.y[:, idx]
            
            if len(idx) == 1:
                if self.yaxis == "abs:deg":
                    out = self.y[:, idx[0]]
                elif self.yaxis == "re:im":
                    out = self.y[:, idx[0]]
            else: 
                if self.yaxis == "abs:deg":
                    out = self.y[:, idx[0]] * np.exp(1j*math.pi*self.y[:, idx[1]]/180.0)
                elif self.yaxis == "re:im":
                    out = self.y[:, idx[0]] + 1j*self.y[:, idx[1]]

            if out == None:
                out = self.y[:, idx]

            if out.size == 1:
                return out[0].squeeze()
            else:
                return out.squeeze()
        else:
            raise  pkex.BasePyKatException("No output by the name '{0}' found in the output".format(str(value)))
      
class katRun2D(object):
    def __init__(self):
        self.runtime
        self.startDateTime = datetime.datetime.now()
        self.x = None
        self.y = None
        self.z = None
        self.xlabel = None
        self.ylabel = None
        self.zlabels = None
        self.katScript = None
        self.katVersion = None
        
    def saveKatRun(self, filename):
        with open(filename,'w') as outfile:
            pickle.dump(self, outfile)
    
    @staticmethod
    def loadKatRun(filename):
        with open(filename,'r') as infile:
            return pickle.load(infile)
    
    def get(self, value): return self[value].squeeze()
    
    def __getitem__(self, value):
        idx = [i for i in range(len(self.zlabels)) if self.zlabels[i].split()[0] == str(value)]
        
        if len(idx) > 0:
            return self.z[idx].squeeze()
        else:
            raise  pkex.BasePyKatException("No output by the name {0} found".format(str(value)))
    
        
class Signals(object):
    class fsig(object):
        def __init__(self, param, name, amplitude, phase, signal):
            self._params = []
            self.__target = param
            self.__name = name
            self.__amplitude = Param("amp", self, SIfloat(amplitude))
            self.__phase = Param("phase", self, SIfloat(phase))
            self.__removed = False
            self.__signal = signal
            
            # unfortunatenly the target names for fsig are not the same as the
            # various parameter names of the components, e.g. mirror xbeta is x 
            # for fsig. So we need to check here what type of component we are targetting
            # and then based on the parameter specfied get the name
            if not param.canFsig:
                raise  pkex.BasePyKatException("Cannot fsig parameter {1} on component {0}".format(str(param._owner().name), param.name))
            
        def _register_param(self, param):
            self._params.append(param)
        
        @property
        def removed(self): return self.__removed
  
        def remove(self):
            self.__signal._kat.remove(self)
            
        def _on_remove(self):
            if self.__removed:
                raise pkex.BasePyKatException("Signal {0} has already been marked as removed".format(self.name))
            else:
                self.__signal.targets.remove(self)
                self.__remove = True
        
        @property
        def name(self): return self.__name

        @property
        def amplitude(self): return self.__amplitude
        @amplitude.setter
        def amplitude(self,value): self.__amplitude.value = SIfloat(value)


        @property
        def phase(self): return self.__phase
        @phase.setter
        def phase(self,value): self.__phase.value = SIfloat(value)

        @property
        def target(self): return self.__target.fsig_name

        @property
        def owner(self): return self.__target._owner().name
    
        def getFinesseText(self):
            rtn = []
    
            for p in self._params:
                rtn.extend(p.getFinesseText())
        
            return rtn
    
    @property
    def name(self):
        # if we don't have any signals yet then use a dummy name
        # however we need to always tune a real fsig command
        # so need to get the name of at least one of them
        # as if you tune one you tune them all
        if len(self.targets) == 0:
            return "fsignal"
        else:
            return self.targets[0].name
            
    @property
    def removed(self): return False # we can never remove the Signal object altogethr just the individual fsig targets

    def remove(self):
        for t in self.targets:
            t.remove()
        
        del self.targets[:]
        
    @property
    def f(self): return self.__f
    @f.setter
    def f(self,value): self.__f.value = SIfloat(value)
    
    def __init__(self, kat):
        self.targets = []
        self._params = []
        self.__f = Param("f", self, 1)
        self._kat = kat
        
    def _register_param(self, param):
        self._params.append(param)
        
    def apply(self, target, amplitude, phase, name=None):
        
        if target == None:
            raise  pkex.BasePyKatException("No target was specified for signal to be applied")
        
        if name == None:
            name = "sig_" + target._owner().name + "_" + target.name
        
        self.targets.append(Signals.fsig(target, name, amplitude, phase, self))
        
    def getFinesseText(self):
        rtn = []
        
        for t in self.targets:
            rtn.extend(t.getFinesseText())
            
            rtn.append("fsig {name} {comp} {target} {frequency} {phase} {amplitude}"
                            .format(name = t.name,
                                    comp=t.owner,
                                    target=t.target,
                                    frequency=str(self.f),
                                    phase=str(t.phase),
                                    amplitude=str(t.amplitude if t.amplitude != None else "")))

        for p in self._params:
            rtn.extend(p.getFinesseText())
        
        return rtn
        
class Block:
    def __init__(self, name):
        self.__name = name
        self.contents = [] # List of objects and strings of finesse code
        self.enabled = True 
        
    @property
    def name(self): return self.__name
    
Constant = namedtuple('Constant', 'name, value, usedBy')

id___ = 0

class kat(object):  

    def __new__(cls, *args, **kwargs):
        # This may seem like an arbitrary step but here we are creating a
        # new class that is a base class of itself. This is because when
        # the kat object adds new components it also adds properties for
        # each of these. There properties are unique to each kat object,
        # but properties are part of the class definition. Thus if two
        # kat objects share the same class definition they also have the
        # same properties regardless of whether they have the actual
        # object added to it. So we create an instance specific class.
        global id___
        id___ += 1
        cnew = type(pykat.finesse.kat.__name__ + str("_") + str(id___), (pykat.finesse.kat,), {})
        return object.__new__(cnew)
    
    def __init__(self, kat_file=None, kat_code=None, katdir="", katname="", tempdir=None, tempname=None):
        self.scene = None # scene object for GUI
        self.verbose = True
        self.__blocks = OrderedDict() # dictionary of blocks that are used
        self.__components = {}  # dictionary of optical components      
        self.__detectors = {}   # dictionary of detectors
        self.__commands = {}    # dictionary of commands
        self.__gui = None
        self.nodes = NodeNetwork(self)  
        self.__katdir = katdir
        self.__katname = katname
        self.__tempdir = tempdir
        self.__tempname = tempname
        self.pykatgui = None
        self.__signals = Signals(self)
        self.constants = {}
        self.vacuum = []
        self.__prevrunfilename = None
        self.printmatrix = None
        
        # initialise default block
        self.__currentTag= NO_BLOCK
        self.__blocks[NO_BLOCK] = Block(NO_BLOCK)
        
        # Various options for running finesse, typicaly the commands with just 1 input
        # and have no name attached to them.
        self.retrace = None
        self.deriv_h = None
        self.scale = None
        self.__trace = None
        self.__phase = None
        self.__maxtem = None
        self.__noxaxis = False
        self.__time_code = None
        self.__yaxis = "abs" # default yaxis
        self.__lambda0 = 1064e-9
        
        if kat_code != None and kat_file != None:
            raise pkex.BasePyKatException("Specify either a Kat file or some Kat code, not both.")
        
        if kat_code != None:
            self.parseCommands(kat_code)
        
        if kat_file != None:
            self.loadKatFile(kat_file)

    def __deepcopy__(self, memo):
        """
        When deep copying a kat object we need to take into account
        the instance specific properties. This is because when
        the kat object adds new components it also adds properties for
        each of these. There properties are unique to each kat object,
        but properties are part of the class definition. Thus if two
        kat objects share the same class definition they also have the
        same properties regardless of whether they have the actual
        object added to it. So we create an instance specific class.
        """
        result = self.__class__.__new__(self.__class__.__base__)
        memo[id(self)] = result
        result.__dict__ = copy.deepcopy(self.__dict__, memo)

        # Find all properties in class we are copying
        # and deep copy these to the new class instance
        for x in self.__class__.__dict__.items():
            if isinstance(x[1], property):
                setattr(result.__class__, x[0], x[1])
    
        result.nodes._NodeNetwork__update_nodes_properties()
                
        # Update any weakrefs
        for c in result.components:
            result.components[c]._Component__update_node_setters()
        
        return result
    
    @property
    def signals(self): return self.__signals

    yaxis_options = ["abs:deg","db:deg","re:im","abs","db","deg"]
    @property
    def yaxis(self): return self.__yaxis
    @yaxis.setter
    def yaxis(self, value):
        
        if not str(value) in self.yaxis_options:
            raise pkex.BasePyKatException("yaxis value '{0}' is not a value option. Valid options are: {1}".format(str(value), ",".join(self.yaxis_options) ))
            
        self.__yaxis = str(value)

    @property
    def trace(self): return self.__trace
    @trace.setter
    def trace(self, value):
        value = int(value)

        if value < 0 or value > 255:
            raise pkex.BasePyKatException('trace command only accepts values in the range 0-255.')
        else:
            self.__trace = value

    @property
    def lambda0(self): return self.__lambda0
    @lambda0.setter
    def lambda0(self, value):
        self.__lambda0 = SIfloat(value)
        
        for node in self.nodes.getNodes():
            if self.nodes[node].q != None:
                self.nodes[node].q.wavelength = self.__lambda0

    @property
    def maxtem(self): return self.__maxtem
    @maxtem.setter
    def maxtem(self,value):
        if value == "off":
            self.__maxtem = -1
        else:
            self.__maxtem = int(value)
    
    @property
    def phase(self): return self.__phase
    @phase.setter
    def phase(self,value): self.__phase = int(value)
        
    @property
    def timeCode(self): return self.__time_code
    @timeCode.setter
    def timeCode(self,value): self.__time_code = bool(value)
    
    @property
    def components(self):
        return self.__components.copy()
    
    @property
    def detectors(self):
        return self.__detectors.copy()
        
    @property
    def noxaxis(self): return self.__noxaxis
    @noxaxis.setter
    def noxaxis(self,value): self.__noxaxis = bool(value) 

    @staticmethod
    def logo():
        print ("""                                              ..-
    PyKat {0:7}         _                  '(
                          \\`.|\\.__...-\"\"""-_." )
       ..+-----.._        /  ' `            .-'
   . '            `:      7/* _/._\\    \\   (
  (        '::;;+;;:      `-"' =" /,`"" `) /
  L.        \\`:::a:f            c_/     n_'
  ..`--...___`.  .    ,  
   `^-....____:   +.      {1}\n""".format(pykat.__version__, pykat_web))
    
    def loadKatFile(self, katfile, blocks=None):
        commands=open(katfile).read()
        self.parseCommands(commands, blocks=blocks)
    
    def parseKatCode(self, code, blocks=None):
        #commands = code.split("\n")
        self.parseCommands(code, blocks=blocks)

    def processConstants(self, commands):
        """
        Before fully parsing a bunch of commands firstly any constants or variables
        to be recorded and replaced.
        """
        
        try:
            constants = self.constants
        
            for line in commands:
                values = line.split()
            
                if len(values)>0 and values[0] == 'const':
                
                    if len(values) >= 3:
                        if values[1] in constants:
                            raise pkex.BasePyKatException('const command with the name "{0}" already used'.format(values[1]))
                        else:
                            constants[str(values[1])] = Constant(values[1], values[2], [])
                    else:
                        raise pkex.BasePyKatException('const command "{0}" was not the correct format'.format(line))
        
            commands_new = []
        
            for line in commands:
                values = line.split()
            
                if len(values) > 0 and values[0] != 'const':
                    # check if we have a var/constant in this line
                    if line.find('$') >= 0:
                        for key in constants.keys():
                            # TODO: need to fix this for checking mulitple instances of const in a single line
                        
                            chars = [' ', '+', '-', '*', '/', ')']

                            for c in chars:
                                none_found = False
                            
                                while not none_found:
                                    if line.find('$'+key+c) > -1:
                                        constants[key].usedBy.append(line)
                                        line = line.replace('$'+key+c, str(constants[key].value)+ c)
                                    else:
                                        none_found = True
                        
                            if line.endswith('$'+key):
                                constants[key].usedBy.append(line)
                                line = line.replace('$'+key, str(constants[key].value))
                        
                    commands_new.append(line)
    
            self.constants = constants
        
            return commands_new
            
        except pkex.BasePyKatException as ex:
            pkex.PrintError("Error processing constants:", ex)
            sys.exit(1)
        
    def parseCommands(self, commands, blocks=None):
        try:
            blockComment = False
        
            commands=self.remove_comments(commands)
        
            commands=self.processConstants(commands)
        
            after_process = [] # list of commands that should be processed after 
                               # objects have been set and created
        
            for line in commands:
                if len(line.strip()) >= 2:
                    line = line.strip()

                    # Looking for block start or end
                    values = line.split()
            
                    if values[0] == "%%%":
                        if values[1] == "FTblock":
                            newTag = values[2]
                    
                            if self.__currentTag != None and self.__currentTag != NO_BLOCK: 
                                warnings.warn("found block {0} before block {1} ended".format(newTag, self.__currentTag))    
                        
                            if newTag in self.__blocks:
                                raise pkex.BasePyKatException("Block `{0}` has already been read".format(newTag))
                        
                            self.__blocks[newTag] = Block(newTag) # create new list to store all references to components in block
                            self.__currentTag = newTag
                    
                        if values[1] == "FTend":
                            self.__currentTag = NO_BLOCK
                    
                        continue

                    # only include listed blocks, if we have specfied them
                    if blocks != None and self.__currentTag not in blocks:
                        continue
                
                    # don't read comment lines
                    if line[0] == "#" or line[0] == "%":
                        continue
            
                    # check if block comment is being used
                    if not blockComment and line[0:2] == "/*":
                        blockComment = True
                        continue
                    elif blockComment and line[0:2] == "*/":
                        blockComment = False
                        continue
            
                    first = line.split(" ",1)[0]
                    obj = None
                
                    if(first == "m" or first == "m1" or first == "m2"):
                        obj = pykat.components.mirror.parseFinesseText(line)
                    elif(first == "s"):
                        obj = pykat.components.space.parseFinesseText(line)
                    elif(first == "l"):
                        obj = pykat.components.laser.parseFinesseText(line)
                    elif(first == "sq"):
                        obj = pykat.components.squeezer.parseFinesseText(line)
                    elif(first[0:2] == "bs"):
                        obj = pykat.components.beamSplitter.parseFinesseText(line)
                    elif(first[0:2] == "gr"):
                        obj = pykat.components.grating.parseFinesseText(line)
                    elif(first[0:4] == "isol"):
                        obj = pykat.components.isolator.parseFinesseText(line)
                    elif(first[0:4] == "lens"):
                        obj = pykat.components.lens.parseFinesseText(line)
                    elif(first[0:3] == "mod"):
                        obj = pykat.components.modulator.parseFinesseText(line)
                    elif(first[0:2] == "ad"):
                        obj = pykat.detectors.ad.parseFinesseText(line)
                    elif(first[0:2] == "bp"):
                        obj = pykat.detectors.bp.parseFinesseText(line)
                    elif(first[0:4] == "gouy"):
                        obj = pykat.detectors.gouy.parseFinesseText(line)
                    elif(first[0:2] == "pd" and first != "pdtype"):
                        obj = pykat.detectors.pd.parseFinesseText(line)
                    elif(first == "qshot" or first == "qshotS" or first == "qshotN"):
                        obj = pykat.detectors.qshot.parseFinesseText(line)
                    elif(first == "qnoised" or first == "qnoisedS" or first == "qnoisedN"):
                        obj = pykat.detectors.qnoised.parseFinesseText(line)
                    elif(first == "xaxis" or first == "xaxis*"):
                        obj = pykat.commands.xaxis.parseFinesseText(line)
                    elif(first[0:2] == "hd"):
                        obj = pykat.detectors.hd.parseFinesseText(line)
                    elif(first.startswith("qhd")):
                        obj = pykat.detectors.qhd.parseFinesseText(line)
                    elif(first == "x2axis" or first == "x2axis*"):
                        obj = pykat.commands.x2axis.parseFinesseText(line)
                    elif(first == "gauss" or first == "gauss*" or first == "gauss**"):
                        after_process.append(line)
                    elif(first == "scale"):
                        after_process.append(line)
                    elif(first == "pdtype"):
                        after_process.append(line)
                    elif(first == "attr"):
                        after_process.append(line)
                    elif(first == "noxaxis"):
                        self.noxaxis = True
                    elif(first == "lambda"):
                        v = line.split()
                        self.lambda0 = SIfloat(v[-1])
                    elif(first == "yaxis"):
                        v = line.split()
                
                        self.yaxis = v[-1]
                    elif(first == "phase"):
                        v = line.split()
                        if len(v) != 2:
                            raise pkex.BasePyKatException("phase command `{0}` is incorrect.".format(line))
                        else:
                            self.phase = int(v[1])
                    elif(first == "maxtem"):
                        v = line.split()
                        if len(v) != 2:
                            raise pkex.BasePyKatException("maxtem command `{0}` is incorrect.".format(line))
                        else:
                            if v[1] == "off":
                                self.maxtem = -1
                            else:
                                self.maxtem = int(v[1])
                    elif(first == "trace"):
                        v = line.split()
                        if len(v) > 2:
                            raise pkex.BasePyKatException("Trace command `{0}` is incorrect.".format(line))
                        elif len(v) == 2:
                            self.trace = v[1]
                    elif(first == "retrace"):
                        v = line.split()
                        if len(v) > 2:
                            raise pkex.BasePyKatException("Retrace command `{0}` is incorrect.".format(line))
                        elif len(v) == 2:
                            self.retrace = v[1]                        
                    elif(first == "deriv_h"):
                        v = line.split()
                        if len(v) != 2:
                            raise pkex.BasePyKatException("deriv_h command `{0}` is incorrect.".format(line))
                        else:
                            self.deriv_h = float(v[1])
                    elif(first == "gnuterm" or first == "pyterm"):
                        if self.verbose:
                            print ("Ignoring Gnuplot/Python terminal command '{0}'".format(line))
                    elif(first == "fsig"):
                        after_process.append(line)
                    elif(first == "noplot"):
                        obj = line
                        self.__blocks[self.__currentTag].contents.append(line) 
                    else:
                        if self.verbose:
                            print ("Parsing `{0}` into pykat object not implemented yet, added as extra line.".format(line))
                    
                        obj = line
                        # manually add the line to the block contents
                        self.__blocks[self.__currentTag].contents.append(line) 
            
                    if obj != None and not isinstance(obj, six.string_types):
                        if self.hasNamedObject(obj.name):
                            getattr(self, obj.name).remove()
                            print ("Removed existing object '{0}' of type {1} to add line '{2}'".format(obj.name, obj.__class__, line))
                    
                        self.add(obj)
                
                
            # now process all the varous gauss/attr etc. commands which require
            # components to exist first before they can be processed
            for line in after_process:
                first = line.split(" ",1)[0]            
                if first == "gauss" or first == "gauss*" or first == "gauss**":
                    pykat.commands.gauss.parseFinesseText(line, self)
                elif (first == "scale"):
                    v = line.split()
                    accepted = ["psd","psd_hf","asd","asd_hf","meter", "ampere", "degs"]
                
                    if len(v) == 3:
                        component_name = v[2]
                    
                        if v[1].lower() in accepted:
                            val = v[1]
                        else:
                            try:
                                val = SIfloat(v[1])
                            except ValueError as ex:
                                raise pkex.BasePyKatException("Line `{0}`:\nAccepted scale values are decimal numbers or %s." % (line,str(accepted)))
                            
                        if component_name in self.__detectors :
                            self.__detectors[component_name].scale.append(val)
                        else:
                            raise pkex.BasePyKatException("scale command `{0}` refers to non-existing output".format(component_name))
                    elif len(v) == 2:
                        if v[1] == "meter" or v[1] == "ampere" or v[1] == "deg":
                            self.scale = v[1]
                        else:
                            self.scale = SIfloat(v[1])
                    else:
                        raise pkex.BasePyKatException("scale command `{0}` is incorrect.".format(line))
                elif (first == "pdtype"):
                    v = line.split()
                    if len(v) == 3:
                        component_name = v[1]
                        if component_name in self.__detectors :
                            self.__detectors[component_name].pdtype = v[2]
                        else:
                            raise pkex.BasePyKatException("pdtype command `{0}` refers to non-existing detector".format(component_name))
                    else:
                        raise pkex.BasePyKatException("pdtype command `{0}` is incorrect.".format(line))
                elif(first == "attr"):
                    v = line.split()

                    if len(v) < 4:
                        raise pkex.BasePyKatException("attr command `{0}` is incorrect.".format(line))
                    else:
                        # get the component/detector in question
                        if v[1] in self.__components:
                            comp = self.__components[v[1]]
                        elif v[1] in self.__detectors:
                            comp = self.__detectors[v[1]]
                        else:
                            raise pkex.BasePyKatException("Could not find the component '{0}' for attr command in line '{1}'".format(v[1], line))
                
                        if len(v[2:]) % 2 == 1:
                            raise pkex.BasePyKatException("Attr command '{0}' must specify both parameter and value pairs".format(line))
                                                
                        # convert split list to key value pairs
                        #kv = dict(itertools.izip_longest(*[iter(v[2:])] * 2, fillvalue=None))
                        kv = dict(izip_longest(*[iter(v[2:])] * 2, fillvalue=None))

                        comp.parseAttributes(kv)
                    
                elif(first == "fsig"):
                
                    v = line.split()
                
                    name = str(v[1])
                
                    if v[2] not in self.__components:
                        raise pkex.BasePyKatException("Could not find the component '{0}'. Line: '{1}'".format(v[2], line))
                
                    comp = self.__components[v[2]]
                
                    if comp._default_fsig() == None:
                        raise pkex.BasePyKatException("Component '{0}' cannot be fsig'd. Line: '{1}'".format(comp.name, line))
                    
                    param = None
                    amp = None
                
                    if len(v) == 5:
                        param == None
                        freq = float(v[3])
                        phase = float(v[4])
                    elif len(v) == 6:
                        if v[3].isdigit():
                            freq = float(v[3])
                            phase = float(v[4])
                            amp = float(v[5])
                        else:
                            param = v[3]
                            freq = float(v[4])
                            phase = float(v[5])
                    elif len(v) == 7:
                        param = v[3]
                        freq = float(v[4])
                        phase = float(v[5])
                        amp = float(v[6])
                    else:
                        raise pkex.BasePyKatException("'{0}' isnot a valid fsig command".format(line))
                    
                    self.signals.f = freq
                    self.signals.apply(comp._default_fsig(), amp, phase, name)
                
                else:
                    raise pkex.BasePyKatException("Haven't handled parsing of '{0}'".format(line))
                    
            self.__currentTag = NO_BLOCK 
        

        except pkex.BasePyKatException as ex:
            pkex.PrintError("Error parsing line: '%s':"%  line, ex)
            sys.exit(1)
            
    def saveScript(self, filename=None):
        """
        Saves the current kat object to a Finesse input file
        """
        try:
            katScript = "".join(self.generateKatScript())       
            katfile = open(filename,'w')
            katfile.writelines(katScript)
            katfile.flush()
            katfile.close()

        except pkex.BasePyKatException as ex:
            print (ex)

            
    def getProcess(self, callback, **kwargs):
        """
        """
        
        cmd = "\n".join(self.generateKatScript())
        
        return Process(target=f__lkat_process, args=(callback, cmd, kwargs))

    def run(self, printout=0, printerr=0, plot=None, save_output=False, save_kat=False, kat_name=None, cmd_args=None, getTraceData=False):
        """ 
        Runs the current simulation setup that has been built thus far.
        It returns a katRun or katRun2D object which is populated with the various
        data from the simulation run.
        printout=1 prints the Finesse banner
        printerr shows the Finesse progress (set kat.verbose=1 to see warnings and errors)
        plot (string) - Sets gnuterm for plotting
        save_output (bool) - if true does not delete out file
        save_kat (bool) - if true does not delete kat file
        kat_name (string) - name of kat file if needed, will be randomly generated otherwise
        cmd_args (list of strings) - command line flags to pass to FINESSE
        getTraceData (bool) - If true a list of dictionaries is returned along with the
        output file. Each dictionary is the result of the beam tracing
        that Finesse performs, the keys are the node names and the values
        are the x and y beam parameters. If no tracing is done a None
        is returned.
        """
        start = datetime.datetime.now()
        
        try:        
            if not hasattr(self, "xaxis") and self.noxaxis != None and self.noxaxis == False:
                raise pkex.BasePyKatException("No xaxis was defined")
            
            if len(self.__katdir) == 0:
                # Get the environment variable for where Finesse is stored
                self.__finesse_dir = os.environ.get('FINESSE_DIR')
                
                if self.__finesse_dir == None :
                    raise pkex.MissingFinesseEnvVar()
            else:
                self.__finesse_dir = self.__katdir
                
            if len(self.__katname) == 0:
                katexe = "kat"
                
                if os.sys.platform == "win32":
                    katexe += ".exe"
            else:
                katexe = self.__katname
            
            kat_exec = os.path.join(self.__finesse_dir, katexe) 
            
            # check if kat file exists and it is executable by user        
            if not (os.path.isfile(kat_exec) and os.access(kat_exec, os.X_OK)):
                raise pkex.MissingFinesse()
                
            if self.verbose: print ("--------------------------------------------------------------")
            if self.verbose: print ("Running kat - Started at " + str(start))
            
            if hasattr(self, "x2axis") and self.noxaxis == False:
                r = katRun2D()
            else:
                r = katRun()
                
            r.yaxis = self.yaxis
            
            r.katScript = "".join(self.generateKatScript())
            r.katScript += "time\n"

            if (plot==None):
                # ensure we don't do any plotting. That should be handled
                # by user themselves
                r.katScript+=("gnuterm no\n")
                r.katScript+=("pyterm no\n")
            else:
                r.katScript+=(plot+"\n")
            
            # create a kat file which we will write the script into
            if self.__tempname == None:
                katfile = tempfile.NamedTemporaryFile(mode ='w', suffix=".kat", dir=self.__tempdir)
            else:
                filepath =os.path.join(self.__tempdir, self.__tempname+".kat" )
                katfile = open( filepath, 'w' ) 
                
            katfile.writelines(r.katScript)
            #katfile.writelines(bytes(r.katScript, 'UTF-8'))
            katfile.flush()

            if printout == 1 or plot != None:
                cmd=[kat_exec]
            else:
                cmd=[kat_exec, '--perl1']
            
            if self.__time_code:
                cmd.append('--perf-timing')
            
            if cmd_args != None:
                cmd.extend(cmd_args)

            if getTraceData:
                cmd.append('--trace')
                
            cmd.append('--no-backspace')
            # set default format so that less repeated numbers are printed to the
            # output file, should speed up running and parsing of output files
            cmd.append('-format=%.15g')

            cmd.append(katfile.name)
                            
            p=subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            err = ""
            
            #if self.verbose: print "Finesse output:"            
            for aline in iter(p.stderr.readline, b""):
                if six.PY2:
                    line = unicode(aline, "utf-8")
                else:
                    line = aline
                    
                if len(line) > 0:
                    # remove any ANSI commands
                    #ansi = re.compile(r'\x1b[^m]*m')
                    #line = ansi.sub('', line)
                    line = re.sub(br'\x1b[^m]*m', '', line, re.UNICODE)

                    # warnings and errors start with an asterisk 
                    # so if verbose show them
                    if line.lstrip().startswith(b'*PROG*'):
                        line = line[8:-1]
                        vals = line.split(b"-",1)
                        action = vals[0].strip()
                        prc = vals[1].strip()[:]
                    
                        if printerr == 1:
                            if six.PY2:
                                sys.stdout.write("\r{0} {1}".format(action, prc))
                            else:
                                sys.stdout.write("\r{0} {1}".format(str(action, 'utf-8'), str(prc, 'utf-8')))
                                
                    elif line.lstrip().startswith(b'*'):
                        if self.verbose:
                            if six.PY2:
                                sys.stdout.write(line)        
                            else:
                                sys.stdout.write(str(line,'utf-8')) 
                                
                    elif line.rstrip().endswith(b'%'):
                        vals = line.split("-")
                        action = vals[0].strip()
                        prc = vals[1].strip()[:]
                        
                        if printerr == 1:
                            if six.PY2:
                                sys.stdout.write("\r{0} {1}".format(action, prc))
                            else:
                                sys.stdout.write("\r{0} {1}".format(str(action, 'utf-8'), str(prc, 'utf-8')))
                            
                    else:
                        if six.PY2:
                            err="".join((err,line))
                        else:
                            err="".join((err,str(line, 'utf-8')))

            
            [out,errpipe] = p.communicate()

            _out = str(out).split("\n")

            for line in _out[::-1]:
                if line.lstrip().startswith('computation time:'):
                    r.runtime = float(line.split(":")[1].replace("s",""))

            if printout == 1: 
                print (out)
            else:
                if printerr == 1: print ("")

            # get the version number
            ix = out.find(b'build ') + 6
            ix2 = out.find(b')',ix)
            r.katVersion = out[ix:ix2]
            
            r.runDateTime = datetime.datetime.now()

            # If Finesse returned an error, just print that and exit!
            if p.returncode != 0:
                raise pkex.FinesseRunError(err, katfile.name)
            
            self.__prevrunfilename = katfile.name
            
            root = os.path.splitext(katfile.name)
            base = os.path.basename(root[0])
            path = os.path.split(katfile.name)[0]            
            outfile = root[0] + ".out"

            traceData = None
            
            if getTraceData:
                # First see if we have any trace files
                
                traceFiles = [file for file in os.listdir(path) if file.endswith(".trace") and file.startswith(base)]
                
                print("Found %i trace files" % len(traceFiles))
                print(path)
                print(traceFiles)
                
                if len(traceFiles) > 0:
                    import fileinput
                    traceData = []
                    
                    for file in traceFiles:
                        traceData.append({})
                        ifile = fileinput.input(os.path.join(path, file))
                    
                        for line in ifile:
                            line = line.strip()
                        
                            if len(line) > 0:
                                a = line.split(':', 1)
                        
                                if a[0].isdigit():
                                    print("Found %s" % a[0])
                                
                                    values = a[1].split()
                                
                                    node_name = values[1].split("(")[0]
                                
                                    line1x = ifile.readline()
                                    line2x = ifile.readline()
                                    line1y = ifile.readline()
                                    line2y = ifile.readline()

                                    qx = line2x.strip().split()[0].split("=")[1]
                                    qy = line2y.strip().split()[0].split("=")[1]
                                
                                    traceData[-1][node_name] = (pykat.beam_param(q=complex(qx)), pykat.beam_param(q=complex(qy)))

            
            if save_output:        
                newoutfile = "{0}.out".format(base)
                
                cwd = os.path.os.getcwd()
                newoutfile = os.path.join(cwd,newoutfile)
                
                if os.path.isfile(newoutfile):
                    os.remove(newoutfile)
                    
                os.rename(outfile, newoutfile)

                if self.verbose: print ("\nOutput data saved to '{0}'".format(newoutfile))

            if len(self.detectors.keys()) > 0:
                if hasattr(self, "x2axis") and self.noxaxis == False:
                    [r.x,r.y,r.z,hdr] = self.readOutFile(outfile)
                
                    r.xlabel = hdr[0]
                    r.ylabel = hdr[1]
                    r.zlabels = [s.strip() for s in hdr[2:]]
                    #r.zlabels = map(str.strip, hdr[2:])
                else:
                    [r.x,r.y,hdr] = self.readOutFile(outfile)
            
                    r.xlabel = hdr[0]
                    r.ylabels = [s.strip() for s in hdr[1:]]
                    #r.ylabels = map(str.strip, hdr[1:]) // replaced 090415 adf 
                    
            if save_kat:
                if kat_name == None:
                    kat_name = "pykat_output"                
                
                cwd = os.path.os.getcwd()
                newkatfile = os.path.join(cwd, kat_name + ".kat")
                
                if os.path.isfile(newkatfile):
                    os.remove(newkatfile)
                  
                os.rename(katfile.name, newkatfile)         
                
                if self.verbose: print ("Kat file saved to '{0}'".format(newkatfile))
                
            if self.trace != None and self.trace > 0:
                #print ("{0}".format(out))
                #if self.trace & 1:
                    #search = out.find(' --- highest order of TEM modes')
                    #if search > -1:
                        #print ("Trace 1: {0}".format(out[search:]))

                # for now, just try to print the trace block in full
                print (out[out.find(' ---') :])

            katfile.close()
            perfData = []

            rtn = [r]
            
            if self.__time_code:
                perffile = open(root[0] + ".perf",'r')
                
                for l in perffile.readlines():
                    vals = l.strip().split()
                    perfData.append((vals[0], long(vals[1]), long(vals[2])))
                    #perfData.append((vals[0], float(vals[1]), float(vals[2]), float(vals[3])))
                    
                rtn.append(perfData)
            
            if getTraceData:
                rtn.append(traceData)
                
            if len(rtn) == 1:
                return rtn[0]
            else:
                return rtn
        except KeyboardInterrupt as ex:
            print("Keyboard interrupt caught, stopped simulation.")
        except pkex.FinesseRunError as ex:
            pkex.PrintError("Error from Finesse:", ex)
        except pkex.BasePyKatException as ex:
            pkex.PrintError("Error from pykat:", ex)
        finally:
            if self.verbose: print ("")
            if self.verbose: print ("Finished in " + str(datetime.datetime.now()-start))
            
    def remove(self, obj):
        try:
            if not isinstance(obj, pykat.finesse.Signals) and not (obj.name in self.__components  or obj.name in self.__detectors or obj.name in self.__commands or obj in self.signals.targets):
                raise pkex.BasePyKatException("{0} is not currently in the simulation".format(obj.name))
            
            if obj.removed:
                raise pkex.BasePyKatException("{0} has already been removed".format(obj.name))        

            nodes = None
        
            # store nodes that this componet is attached to as a reference for gui
            if isinstance(obj, Component):
                nodes = self.nodes.getComponentNodes(obj)

            if isinstance(obj, Component):    
                del self.__components[obj.name]
                self.__del_component(obj)
                self.nodes._removeComponent(obj)
            elif isinstance(obj, Command):    
                del self.__commands[obj.name]
                self.__del_command(obj)
            elif isinstance(obj, Detector):    
                del self.__detectors[obj.name]
                self.__del_detector(obj)
            elif isinstance(obj, pykat.finesse.Signals):
                obj.remove()
            elif isinstance(obj, pykat.finesse.Signals.fsig):
                obj._on_remove()
            
            for b in self.__blocks:
                if obj in self.__blocks[b].contents:
                    self.__blocks[b].contents.remove(obj)
        
            if self.pykatgui != None:
                self.pykatgui._onComponentRemoved(obj, nodes)
    
            del nodes
        
            #import gc
            #print (gc.get_referrers(obj))
            
        except pkex.BasePyKatException as ex:
            pkex.PrintError("Error on removing object:", ex)
            
    def getMatrices(self):
        
        import scipy
        from scipy.sparse import coo_matrix
        
        prev = self.noxaxis
        
        self.noxaxis = True
        self.printmatrix = True
        print ("".join(self.generateKatScript()))
        self.verbose = True
        self.run(printout=1)
        self.printmatrix = None
        self.noxaxis = prev        
        
        if self.__prevrunfilename == None:
            return None
        else:
            
            Mcarrier = None
            Msignal = None
            
            if os.path.exists("klu_full_matrix_car.dat"):
                M = np.loadtxt("klu_full_matrix_car.dat")
                
                if M.size > 0:
                    row = M[:,0]-1
                    col = M[:,1]-1
                    data = M[:,2] + 1j * M[:,3]
                    N = row.max()+1
                    Mcarrier = coo_matrix((data,(row,col)), shape=(N,N))
                
        
            if os.path.exists("klu_full_matrix_sig.dat"):
                M = np.loadtxt("klu_full_matrix_sig.dat")
                
                if M.size > 0:
                    row = M[:,0]-1
                    col = M[:,1]-1
                    data = M[:,2] + 1j * M[:,3]
                    N = row.max()+1
                    Msignal = coo_matrix((data,(row,col)), shape=(N,N))
        
            return (Mcarrier, Msignal)

    
    def hasNamedObject(self, name):
        return name in self.__components or name in self.__detectors or name in self.__commands
        
    def add(self, obj):
        try:
            obj.tag = self.__currentTag
            self.__blocks[self.__currentTag].contents.append(obj)
            
            if isinstance(obj, Component):
                
                if obj.name in self.__components :
                    raise pkex.BasePyKatException("A component with name '{0}' has already been added".format([obj.name]))            
                            
                self.__components[obj.name] = obj
                self.__add_component(obj)
                
            elif isinstance(obj, Detector):
                
                if obj.name in self.__detectors :
                        raise pkex.BasePyKatException("A detector '{0}' has already been added".format(obj.name))
                        
                self.__detectors[obj.name] = obj
                self.__add_detector(obj)
                
            elif isinstance(obj, Command):
                
                self.__commands[obj.__class__.__name__] = obj
                self.__add_command(obj)
                
            else:
                raise pkex.BasePyKatException("Object {0} could not be added".format(obj))
                
            obj._on_kat_add(self)
            
        except pkex.BasePyKatException as ex:
            pkex.PrintError("Error on adding object:", ex)

    def readOutFile(self, filename):
        
        with open(filename,'r') as outfile:
            # read first to lines to get to header line
            outfile.readline()
            outfile.readline()
            
            hdr = outfile.readline().replace('%','').replace('\n','').split(',')
        
        data = np.loadtxt(filename, comments='%',skiprows=4)

        # convert 1D arrays into 2D ones for simpler selection
        if len(data.shape) == 1:
            data = np.array([data])
                            
        if hasattr(self, "x2axis") and self.noxaxis == False:
            # need to parse 2D outputs slightly different as they are effectively 2D matrices
            # written in linear form
            x = data[0::(1+self.x2axis.steps),0].squeeze()
            y = data[0:(1+self.x2axis.steps),1]
            # get rows and columns lined up so that we can reshape a single column of all x/y data
            # into a matrix
            z = data[:,2:].transpose().reshape(data.shape[1]-2, 1+self.xaxis.steps, 1+self.x2axis.steps).squeeze()
            # once you do this the data for y and x axes need swapping
            z = z.swapaxes(1,2)
            return [x, y, z, hdr]
        else:
            shape_len = len(data.shape)
            
            rows,cols = data.shape
            
            x = data[:,0].squeeze()
            y = data[:,1:cols]
        
            return [x, y, hdr]

    def removeLine(self, fragment) :
        """
        This will search all blocks and search for the string
        fragment specified and remove it.
        WARNING: This will only remove non-parsed commands, it will not
        remove commands that have already been parsed
        into a pykat object, such as mirrors and beamsplitters, use
        kat.remove or kat.component.remove() to delete parsed objects.
        """
        found = False

        for key in self.__blocks:
            objs = self.__blocks[key].contents
            for obj in objs:
                if isinstance(obj,  six.string_types):
                    if fragment in obj:
                        print ("  ** removing line '{0}'".format(obj))
                        objs.remove(obj)
                        found = True

        if not found:
            pkex.BasePyKatException("The command fragment '%s' is not an extra line added to this kat object. Please check that the item you are trying to remove has not been parsed as a pykat object." % fragment)


    def addLine(self, line, block=NO_BLOCK) :
        """
        This will forcefully add a line of FINESSE code to a particular block
        if specfied. This command will not undergo any parsing so it will remain
        as just a string. This of course can create possible conflicts with other
        pykat object that create similar commands so becareful.
        """
        self.__blocks[block].contents.append(line)
        
    def printExtraLines(self):
        """
        This prints all the Finesse commands that have not been parsed
        into pykat objects. This should be used for reference only. To
        add or remove extra lines use the addLine and removeLine methods.
        """
        found = False
        
        for key in self.__blocks:
            objs = self.__blocks[key].contents
            for obj in objs:
                if isinstance(obj, six.string_types):
                    print(obj)
                    found = True
        
        if not found:
            print("No extra lines were found")
        
                        
    def generateKatScript(self) :
        """ Generates the kat file which can then be run """

        def writeBlock():
            for obj in objs:
                if isinstance(obj, six.string_types):
                    out.append(obj + '\n')
                    
                elif isinstance(obj, Component) or isinstance(obj, Detector) or isinstance(obj, Command):
                    txt = obj.getFinesseText() 
                    
                    if txt != None:
                        if isinstance(txt,list):
                            for t in txt:
                                out.append(t + "\n")
                        else:
                            out.append(txt + "\n")

        
        out = []    
        import datetime
        strtoday = datetime.datetime.now()
        out.append(strtoday.strftime("%% Generated by PyKat %d.%m.%Y %H:%M:%S\n") )

        # write the FTblocks
        for key in self.__blocks:
            objs = self.__blocks[key].contents

            if key != NO_BLOCK:
                if np.size(objs)>0:
                    out.append("\n")
                    out.append("%%% FTblock " + key + "\n")
                    writeBlock()
                    out.append("%%% FTend " + key + "\n")

        # write the NO_BLOCK blocks
        for key in self.__blocks:
            objs = self.__blocks[key].contents


            if key == NO_BLOCK:
                if np.size(objs)>0:
                    out.append("\n")
                    writeBlock()
                
        # now loop through all the nodes and get any gauss commands
        for key in self.nodes.getNodes():
            txt = self.nodes.getNodes()[key].getFinesseText()
            
            if txt != None:
                if isinstance(txt,list):
                    for t in txt: out.append(t+ "\n")
                else:
                    out.append(txt + "\n")
        

        # now get any signal commands
        txt = self.signals.getFinesseText()
        
        if txt != None:
            if isinstance(txt,list):
                for t in txt: out.append(t+ "\n")
            else:
                out.append(txt + "\n")

        if self.vacuum != None:
            
            if isinstance(self.vacuum, collections.Container):
                objs = []
                
                if len(self.vacuum) > 0:
                    for a in self.vacuum:
                        if hasattr(a, 'name'):
                            objs.append(a.name)
                        else:
                            objs.append(str(a))

                    out.append("vacuum {0}\n".format(" ".join(objs)))
                                        
            elif isinstance(self.vacuum, six.string_types):
                out.append("vacuum {0}\n".format(self.vacuum))
            else:
                pkex.BasePyKatException("Couldn't understand vacuum input list")

        if self.scale != None and self.scale !='': out.append("scale {0}\n".format(self.scale))
        if self.phase != None: out.append("phase {0}\n".format(self.phase))
        if self.trace != None: out.append("trace {0}\n".format(self.trace))
        if self.maxtem != None:
                if self.maxtem == -1:
                        out.append("maxtem off\n")
                else:
                        out.append("maxtem {0}\n".format(self.maxtem))

        if self.noxaxis == True:
            out.append("noxaxis\n")
            
        if self.yaxis != None:
            out.append("yaxis {0}\n".format(self.yaxis))
         
        if self.printmatrix != None and self.printmatrix == True:
            out.append("printmatrix\n")
            
        if self.lambda0 != 1064e-9:
            out.append("lambda {0}\n".format(self.lambda0))
            
        # ensure we don't do any plotting. That should be handled
        # by user themselves
        #out.append("gnuterm no\n")
        #out.append("pyterm no\n")
        
        return out

    def optivis(self):
        if not HAS_OPTIVIS:
            print("Optivis is not installed")
            return None
        
        import optivis.scene as scene
        import optivis.bench.links as links
        import optivis.view.canvas as canvas
        
        scene = scene.Scene(title="My pykat layout")
        
        # Run through again to add links
        for c in self.getComponents():
            if not isinstance(c, pykat.components.space):
                continue
            
            a = c.connectingComponents()
            
            # Need to handle where spaces don't connect two components but there is a loose
            # node, which may or may not have detectors attached
            if a[0] is None or a[1] is None:
                continue
                
            c1 = a[0].getOptivisComponent()
            c2 = a[1].getOptivisComponent()
            
            no = a[0].getOptivisNode("Output", c.nodes[0])
            ni = a[1].getOptivisNode("Input", c.nodes[1])
            
            if no is None or ni is None:
                raise pkex.BasePyKatException("Optivis node is None")
            
            c._optivis_component = links.Link(outputNode=no, inputNode=ni, length=c.L.value)
            
            c.label_node1 = optivis_label(text="", position=optivis_coord(-0.5, 0), item=c._optivis_component)
            c.label_node2 = optivis_label(text="", position=optivis_coord( 0.5, 0), item=c._optivis_component)
            label_name = optivis_label(text="", position=optivis_coord(0, -0.5), item=c._optivis_component)
            label_name.content["Name"] = c.name
            
            c._optivis_component.labels.append(c.label_node1)
            c._optivis_component.labels.append(c.label_node2)
            c._optivis_component.labels.append(label_name)
            scene.addLink(c._optivis_component)
        
        gui = canvas.Full(scene=scene)
        
        ### get menu bar from Optivis
        menubar = gui.qMainWindow.menuBar()
        
        ### add new Pykat menu and menu items
        pykatMenu = menubar.addMenu('&Pykat')
    
        # save
        save = PyQt4.QtGui.QAction('Save', gui.qMainWindow)
        save.setShortcut('Ctrl+S')
        save.triggered.connect(lambda: self._optivis_doSave(gui))
        pykatMenu.addAction(save)
    
        # trace
        trace = PyQt4.QtGui.QAction('Trace', gui.qMainWindow)
        trace.setShortcut('Ctrl+T')
        trace.triggered.connect(lambda: self._optivis_doTrace(gui))
        pykatMenu.addAction(trace)
    
        return gui
    
    def _optivis_doTrace(self, gui, **kwargs):
        """
        Change at some point to use a stored GUI reference
        """
        if not HAS_OPTIVIS:
            print("Optivis is not installed")
            return None
        
        prev = self.noxaxis
        
        self.noxaxis = True
        out, tdata = self.run(getTraceData=True, **kwargs)
        self.noxaxis = prev
        
        # For now just select the first trace computed
        # Later we could add some gui list to show the different ones
        tdata = tdata[0]
        
        for c in self.getComponents():
            if not isinstance(c, pykat.components.space):
                continue
                
            if not (hasattr(c, "label_node1") and hasattr(c, "label_node2")):
                continue
                
            c.label_node1.content["w0_x"] = tdata[c.nodes[0].name][0].w0
            c.label_node1.content["w_x"] = tdata[c.nodes[0].name][0].w
            c.label_node1.content["z_x"] = tdata[c.nodes[0].name][0].z
            c.label_node1.content["Rc_x"] = tdata[c.nodes[0].name][0].Rc
            c.label_node1.content["Zr_x"] = tdata[c.nodes[0].name][0].zr
            
            c.label_node1.content["w0_y"] = tdata[c.nodes[0].name][1].w0
            c.label_node1.content["w_y"] = tdata[c.nodes[0].name][1].w
            c.label_node1.content["z_y"] = tdata[c.nodes[0].name][1].z
            c.label_node1.content["Rc_y"] = tdata[c.nodes[0].name][1].Rc
            c.label_node1.content["Zr_y"] = tdata[c.nodes[0].name][1].zr
            
            c.label_node2.content["w0_x"] = tdata[c.nodes[1].name][0].w0
            c.label_node2.content["w_x"] = tdata[c.nodes[1].name][0].w
            c.label_node2.content["z_x"] = tdata[c.nodes[1].name][0].z
            c.label_node2.content["Rc_x"] = tdata[c.nodes[1].name][0].Rc
            c.label_node2.content["Zr_x"] = tdata[c.nodes[1].name][0].zr
            
            c.label_node2.content["w0_y"] = tdata[c.nodes[1].name][1].w0
            c.label_node2.content["w_y"] = tdata[c.nodes[1].name][1].w
            c.label_node2.content["z_y"] = tdata[c.nodes[1].name][1].z
            c.label_node2.content["Rc_y"] = tdata[c.nodes[1].name][1].Rc
            c.label_node2.content["Zr_y"] = tdata[c.nodes[1].name][1].zr
       
        gui.redraw()
        
    def _optivis_doSave(self, gui, **kwargs):
        """
        Save kat script from Optivis
        """
        if not HAS_OPTIVIS:
            print("Optivis is not installed")
            return None
    
        # generate file path
        directory = os.path.join(os.path.expanduser('~'), '{0}.kat'.format(gui.scene.title))

        # desired save path
        path = None
    
        # get path to file to export to
        while True:    
            dialog = PyQt4.Qt.QFileDialog(parent=gui.qMainWindow, caption='Save Finesse Script', directory=directory)
            dialog.setAcceptMode(PyQt4.Qt.QFileDialog.AcceptSave)
            dialog.setFileMode(PyQt4.Qt.QFileDialog.AnyFile)

            # show dialog
            dialog.exec_()
      
            if len(dialog.selectedFiles()) is 0:
                # no filename specified
                return

            # get file path and format
            path = str(dialog.selectedFiles()[0])
      
            try:
                # check if we can write to the path
                open(path, 'w').close()
                os.unlink(path)
        
                # if we get this far, all is good so we can break out the infinite loop
                break
            except OSError:
                PyQt4.Qt.QMessageBox.critical(gui.qMainWindow, 'Filename invalid', 'The specified filename is invalid')
            except IOError:
                PyQt4.Qt.QMessageBox.critical(gui.qMainWindow, 'Permission denied', 'You do not have permission to save the file to the specified location.')
    
        # save kat file
        self.saveScript(filename=path)
    
    def openGUI(self):
        if not USE_GUI:
            raise NoGUIException
        else:
            self.app = QCoreApplication.instance() 
            created = False
            
            if self.app == None:
                created = True
                self.app = QApplication([""])
                
            if self.pykatgui == None:
                self.pykatgui = pyKatGUI(self)
                self.pykatgui.main()
            else:
                self.pykatgui.show()
                
            if created: self.app.exec_()
    
    def getComponents(self):
        return self.__components.values()
    
    def hasComponent(self, name):
        return (name in self.__components)
    
    def _newName(self, container, prefix):
        n = 1
        name = "{0}{1}".format(prefix, n)
        
        while name in container:
            n += 1
            name = "{0}{1}".format(prefix,n)
        
        return name
    
    def getNewComponentName(self,prefix):
        '''
        Returns a name for a component which hasn't already been added.
        Returns [prefix] + number, where number is greater than 1. e.g.
        if m1 exists getNewName('m') will return 'm2'
        '''
        return self._newName(self.__components, prefix)
    
    def getNewDetectorName(self,prefix):
        '''
        Returns a name for a component which hasn't already been added.
        Returns [prefix] + number, where number is greater than 1. e.g.
        if m1 exists getNewName('m') will return 'm2'
        '''
        return self._newName(self.__detectors, prefix)
        
    def getNewNodeNames(self,prefix,N=1):
        '''
        Returns a list of names for N number of nodes which haven't already been added.
        Returns [prefix] + number, where number is greater than 1. e.g.
        if m1 exists getNewName('m') will return 'm2'
        '''
        rtn = []
        n = 1
        
        for M in range(1,N+1):
            name = "{0}{1}".format(prefix, n)
            
            while name in self.nodes.getNodes() or (name in rtn):
                n += 1
                name = "{0}{1}".format(prefix,n)
        
            rtn.append(name)
            
        return rtn
        
    def lkat_trace(self, getCavities=True, getNodes=True, getSpaces=True):
        """
        Given the current state of the kat object a new FINESSE process is called and just
        the beam tracing routine is run. The object that is returned contains all the information
        from the beam tracing routine for each node and space components defined as well as cavity 
        commands.   
        """
        if lkat_location == None:
            raise RuntimeError("Could not find shared library 'libkat', please install to a system location or copy to the same directory as this script")
            
        trace_info = Manager().dict()
        
        prev = self.maxtem
        self.maxtem = 0
        
        try:
            p = self.getProcess(f__lkat_trace_callback, trace_info=trace_info,
                                getCavities=getCavities, getNodes=getNodes, getSpaces=getSpaces)
            p.start()
            p.join()
            p.terminate()
            
            # return a local copy of the trace information dictionary
            return dict(trace_info)
        finally: 
            self.maxtem = prev
    
    def __add_detector(self, det):

        if not isinstance(det, Detector):
            raise exceptions.ValueError("Argument is not of type Detector")
        
        name = det.name
        fget = lambda self: self.__get_detector(name)
        
        setattr(self.__class__, name, property(fget))
        setattr(self, '__det_' + name, det)                

    def __del_detector(self, det):

        if not isinstance(det, Detector):
            raise pkex.BasePyKatException("Argument is not of type Detector")
        
        name = det.name
        
        delattr(self.__class__, name)
        delattr(self, '__det_' + name) 
        
    def __get_detector(self, name):
        return getattr(self, '__det_' + name) 
        
    def __add_command(self, com):

        if not isinstance(com, Command):
            raise pkex.BasePyKatException("Argument is not of type Command")
        
        name = com.__class__.__name__
        fget = lambda self: self.__get_command(name)
        
        setattr(self.__class__, name, property(fget))
        setattr(self, '__com_' + name, com)                   

    def __del_command(self, com):

        if not isinstance(com, Command):
            raise exceptions.ValueError("Argument is not of type Command")
        
        name = com.__class__.__name__
        
        #print (getattr(self.__class__, name))
        
        delattr(self.__class__, name)
        delattr(self, '__com_' + name)
        
    def __get_command(self, name):
        return getattr(self, '__com_' + name)            
    
    def __add_component(self, comp):

        if not isinstance(comp, Component):
            raise pkex.BasePyKatException("Argument is not of type Component")
            
        fget = lambda self: self.__get_component(comp.name)
        
        setattr(self.__class__, comp.name, property(fget))
        setattr(self, '__comp_' + comp.name, comp)                   
        
    def __del_component(self, comp):

        if not isinstance(comp, Component):
            raise pkex.BasePyKatException("Argument is not of type Component")
        
        delattr(self.__class__, comp.name)
        delattr(self, '__comp_' + comp.name)
        
    def __get_component(self, name):
        return getattr(self, '__comp_' + name)        

    def remove_comments(self, string):
        """
        This takes a raw Finesse code string and removes any comments
        It returns a list of lines however, not a multiline string.
        Also removes any extrawhite space in command lines.
        """
        pattern = r"(\".*?\"|\'.*?\'|%{3}[^\r\n]*$)|(/\*.*?\*/|%[^\r\n]*$|#[^\r\n]*$|//[^\r\n]*$)"
        # first group captures quoted strings (double or single)
        # second group captures comments (//single-line or /* multi-line */)
        regex = re.compile(pattern, re.MULTILINE|re.DOTALL)
        def _replacer(match):
            # if the 2nd group (capturing comments) is not None,
            # it means we have captured a non-quoted (real) comment string.
            if match.group(2) is not None:
                return "" # so we will return empty to remove the comment
            else: # otherwise, we will return the 1st group
                return match.group(1) # captured quoted-string
        
        # remove any inline comments
        string = regex.sub(_replacer, string)
        
        commands = []
        
        for line in string.split('\n'):
            line = line.replace('\r','')
            if len(line) > 0:
                # remove any mutliple whitespace
                line = " ".join(line.split())
                # add to a list all the positions of any inline comment markers
                i = [line.find('#'), line.find('\\')]
                #i = filter(lambda a: a != -1, i)
                i = [a for a in i if a != -1]    

                if len(i) == 0:
                    commands.append(line)
                else:
                    line = line[0:min(i)]
                    if len(line):
                        commands.append(line)
        
        
        return commands

# printing pykat logo on first input
kat.logo()
