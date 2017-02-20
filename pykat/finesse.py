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

import warnings
import codecs
import uuid
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

from subprocess import Popen, PIPE

try:
    # Python 2
    from itertools import izip_longest
except ImportError:
    # Python 3
    from itertools import zip_longest as izip_longest


try:
    # Add exception in Python 2
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError
    
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
from pykat.external import progressbar
from pykat.freeze import canFreeze
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
    
    if lkat_location is None:
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
                                                 


class BlockedKatFile(object):
    """
    Allows manipulation of blocked kat file.
    
    Example:
        bkf = BlockedKatFile()
        
        bkf.read(katfile)
        bkf.add('tester', "blah\nblah", addAfter="Tunings")
        bkf.remove("Laser")
        bkf.write("mytest.kat")
    """
    
    def __str__(self):
         rtn = ""
         
         for block in self.ordering:
             rtn += "\n%%% FTblock " + block + "\n"
             rtn += self.blocks[block]
             rtn += "%%% FTend " + block + "\n"
         
         return rtn
         
    def __init__(self, NO_BLOCK="NO_BLOCK"):
        self.__NO_BLOCK = NO_BLOCK
        self.ordering = [self.__NO_BLOCK]
        self.blocks = {self.__NO_BLOCK:""}
        self.__currentBlock = self.__NO_BLOCK
        
    def remove(self, *blocks):
        if len(blocks[0]) > 1 and not isinstance(blocks[0], six.string_types):
            # if we've got an iterable thing that isn't a string, eg list or tuple
            # just use that
            blocks = blocks[0]
        
        for block in blocks:
            if block not in self.ordering or block not in self.blocks:
               raise Exception("%s block not found")

            self.ordering.remove(block)
            self.blocks.pop(block)
            
    def add(self, block, contents, addAfter=None):

        if block in self.ordering or block in self.blocks:
            raise Exception("%s block already present")
    
        if addAfter is not None:
            self.ordering.insert(self.ordering.index(addAfter)+1, block)
        else:
            self.ordering.append(block)
            
        self.blocks[block] = contents + "\n"
        
    def write(self, katfile):
        with open(katfile, "w") as f:
            for block in self.ordering:
                f.write("\n%%% FTblock " + block + "\n")
                f.write(self.blocks[block])
                f.write("%%% FTend " + block + "\n")
    
    def read(self, katfile):
        """
        For a given kat file, the blocks are parsed into dictionary as raw strings.
        """
    
        with open(katfile, "r") as f:
            commands = f.readlines()
    
        for line in commands:
            line = line.strip()

            # Looking for block start or end
            values = line.split()
    
            if len(values) >= 3 and values[0] == "%%%":
                if values[1] == "FTblock":
                    newTag = values[2]
                    
                    if self.__currentBlock != None and self.__currentBlock != self.__NO_BLOCK: 
                        warnings.warn("found block {0} before block {1} ended".format(newTag, self.__currentBlock))    
    
                    if newTag in self.blocks:
                        #raise pkex.BasePyKatException("Block `{0}` has already been read".format(newTag))
                        self.__currentBlock = newTag
                        continue
    
                    self.blocks[newTag] = ""
                    self.__currentBlock = newTag
                    self.ordering.append(newTag)
                
                if values[1] == "FTend":
                    self.__currentBlock = self.__NO_BLOCK
        
                continue
                
            if(len(line) == 0 and (self.__currentBlock == self.__NO_BLOCK)):
                continue
        
            self.blocks[self.__currentBlock] += line + "\n"
    
    
class KatBatch(object):
    
    def __init__(self):
        from IPython.parallel import Client

        self._c = Client()
        self._lb = c.load_balanced_view()
        self.lb.block = False
        
        self._todo = []
    
    def _run(dir, commands, **kwargs):
        import pykat
        kat = pykat.finesse.kat()
        kat.verbose = False
        kat.parseCommands(commands)
        
        kw = dict()
        
        if "cmd_args" in kwargs:
            kw["cmd_args"] = kwargs["cmd_args"]
        
        return kat.run(**kw)
    
    def addKat(self, kat, **kwargs):
        import os
        cdir = os.getcwd()
        script = "\n".join(kat.generateKatScript())
        self.todo.append(self.lb.apply_async(self._run, script, **kwargs))
        return self.todo[-1]
    
    def wait(self):
        return self.lb.wait(self.todo)
        
    def results(self):
        return self.todo

                                  
def GUILength(L):
    """
    Should scale the lengths in some way to handle km and mm for time being
    """
    return L # * ( 40 * erfc(L/400.0) + 0.01)

@canFreeze
class KatRun(object):
    def __init__(self):
        self._unfreeze()
        self.runtime = None
        self.StartDateTime = datetime.datetime.now()
        self.x = None
        self.stdout = None
        self.stderr = None
        self.runDateTime = None
        self.y = None
        self.xlabel = None
        self.ylabels = None
        self.katScript = None
        self.katVersion = None
        self.yaxis = None
        self._freeze()
        
    def info(self):
        
        kat = pykat.finesse.kat()
        kat.verbose = False
        kat.parseCommands(self.katScript)
        
        detectors = list(set([lbl.split()[0] for lbl in self.ylabels]))
        detectors.sort()
        
        print("")
        print("--- Output info ---")
        print("")
        print("Run date and time: %s" % self.StartDateTime)
        print("Detectors used: %s" % (", ".join(detectors)))
        print("")

        if kat.noxaxis:
            print("No xaxis used")
        else:
            print("One xaxis used: %s" % kat.xaxis.getFinesseText())
            
        import numpy as np

        maxs = np.max(self.y, 0)
        mins = np.min(self.y, 0)
        
        maxlbl = max([len(lbl) for lbl in self.ylabels])    
        
        for i, lbl in enumerate(self.ylabels):
            a = "{0:" + str(maxlbl) + "} : min = {1:.15e} max = {2:.15e}"
            print(a.format(lbl, mins[i], maxs[i]))
        
        
    def plot(self, detectors=None, filename=None, show=True,
                   yaxis=None, legend=True, loc=0, title=None, styles=None,
                   ylabel=None, y2label=None, xlabel=None, x2label=None,
                   xlim=None, x2lim=None, ylim=None, y2lim=None):
        """
        This will generate a plot for the output data of this particular pykat run.
        It will attempt to generate a plot that shows all the various traces and plots
        by default for quick viewing of the data. Similar to that which would be
        generated by running the Finesse file from the command line.
        
        There are some additional keyword options to customise the plot output slightly:
        
            detectors:          a list of detectors that you want to plot
            filename:           providing a filename here will save the plot to a file.
                                The format is given by the extension provided.
            show:               True | False - whether to display the plot or not
            yaxis:              Set the Finesse yaxis command to base the plot on. By
                                default the original one will be used.
            legend:             True | False - whether to include a legend
            loc:                Location value for the legend, the usual matplotlib one.
            title:              Provide a title for the plot if required.
            styles:             A dictionary which keys being the detector names and the
                                value being a colour and linestyle of the sort 'k:'
            ylabel, xlabel:     Text for the first plot x and y labels
            y2label, x2label:   Text for the second plot x and y labels

            xlim, ylim:         Limits of x- and y-axes of the first plot. List or tuple
                                of length 2.
            x2lim, y2lim:       Limits of x- and y-axes of the second plot. List or tuple
                                of length 2.
        """
        import matplotlib.pyplot as pyplot
        import pykat.plotting as plt
            
        if not show:
            pyplot.ioff()

        kat = pykat.finesse.kat()
        kat.verbose = False
        kat.parseCommands(self.katScript)

        if kat.noxaxis == True:
            raise  pkex.BasePyKatException("This kat object has noxaxis=True, so there is nothing to plot.")

        original_yaxis = kat.yaxis

        if yaxis is not None:
            kat.yaxis = yaxis

        if "log" in kat.yaxis:
            if kat.xaxis.scale == "log":
                plot_cmd = pyplot.loglog
            else:
                plot_cmd = pyplot.semilogy
        else:
            if kat.xaxis.scale == "log":
                plot_cmd = pyplot.semilogx
            else:
                plot_cmd = pyplot.plot

        dual_plot = False
        _func1 = np.abs
        _func2 = None

        plot_cmd1 = None
        plot_cmd2 = None
        
        if "re:im" in kat.yaxis:
            _func1 = np.real
            _func2 = np.imag
            plot_cmd1 = plot_cmd2 = plot_cmd
            
            dual_plot = True
        elif "abs:deg" in kat.yaxis:
            _func1 = np.abs
            _func2 = lambda x: np.rad2deg(np.angle(x))

            plot_cmd1 = plot_cmd
            plot_cmd2 = pyplot.plot if kat.xaxis.scale == "lin" else pyplot.semilogx
            
            dual_plot = True
        elif "db:deg" in kat.yaxis:
            if "db" not in original_yaxis:
                _func1 = lambda x: 10*np.log10(x)
            else:
                _func1 = lambda x: x
                
            _func2 = lambda x: np.rad2deg(np.angle(x))

            plot_cmd1 = plot_cmd
            plot_cmd2 = pyplot.plot if kat.xaxis.scale == "lin" else pyplot.semilogx
            
            dual_plot = True
        elif "abs" in kat.yaxis:
            # _func1 = np.abs
            _func1 = np.real
            plot_cmd1 = plot_cmd
        elif "db" in kat.yaxis:
            if "db" not in original_yaxis:
                _func1 = lambda x: 10*np.log10(x)
            else:
                _func1 = lambda x: x
                
            plot_cmd1 = plot_cmd
        elif "deg" in kat.yaxis:
            _func1 = lambda x: np.rad2deg(np.angle(x))
            plot_cmd1 = plot_cmd
            
        if dual_plot:
            fig = plt.figure(width="full", height=1)
        else:
            fig = plt.figure(width="full")   

        if detectors is None:
            detectors = [lbl.split()[0] for lbl in self.ylabels]
            
        detectors = list(set(detectors))
        detectors.sort()
        
        for det in detectors:
            if not hasattr(kat, det) or (hasattr(kat, det) and not getattr(kat, det).noplot):
                
                if dual_plot:
                    ax = pyplot.subplot(2,1,1)
                    
                if styles is not None and det in styles:
                    l, = plot_cmd1(self.x, _func1(self[det]), styles[det], label=det)
                else:
                    l, = plot_cmd1(self.x, _func1(self[det]), label=det)
                
                if dual_plot: 
                    pyplot.subplot(2,1,2)
                    plot_cmd2(self.x, _func2(self[det]), color=l.get_color(), ls=l.get_linestyle(), label=det)

        if dual_plot:
            if ylabel is None:
                if "abs" in kat.yaxis: ylabel = "Absolute [au]"
                if "re" in kat.yaxis:  ylabel = "Real part [au]"
                    
            if y2label is None:
                if "deg" in kat.yaxis: y2label = "Phase [deg]"
                if "im" in kat.yaxis:  y2label = "Imaginary part [au]"
                    
            if xlim is None:
                xlim = (self.x.min(), self.x.max())
                
            if x2lim is None:
                 x2lim = (self.x.min(), self.x.max())
                 
        else:
            if ylabel is None:
                ylabel = "[au]"
                
            if xlim is None:
                xlim = (self.x.min(), self.x.max())
            
        
        if xlabel is None:
            xlabel = self.xlabel
            
        if x2label is None:
            x2label = self.xlabel
            
        font_label_size = pyplot.rcParams["font.size"]-1
        
        if dual_plot:
            ax = pyplot.subplot(2,1,1)
            pyplot.xlabel(xlabel, fontsize=font_label_size)
            pyplot.ylabel(ylabel, fontsize=font_label_size)
            pyplot.xlim(xlim[0], xlim[1])
            if ylim is not None:
                pyplot.ylim(ylim[0],ylim[1])

            if title is not None:
                pyplot.title(title, fontsize=font_label_size)
    
            pyplot.subplot(2,1,2)
            pyplot.xlabel(x2label, fontsize=font_label_size)
            pyplot.ylabel(y2label, fontsize=font_label_size)
        
            pyplot.xlim(x2lim[0], x2lim[1])
            if y2lim is not None:
                pyplot.ylim(y2lim[0],y2lim[1])
            
        else:
            pyplot.xlabel(xlabel, fontsize=font_label_size)
            pyplot.ylabel(ylabel)
            pyplot.xlim(self.x.min(), self.x.max())
            
            if title is not None:
                pyplot.title(title, fontsize=font_label_size)
            if ylim is not None:
                pyplot.ylim(ylim[0],ylim[1])
    
        pyplot.margins(0, 0.05)
        pyplot.tight_layout()
    
        if legend:
            fig.axes[0].legend(loc=loc, fontsize=font_label_size)
        
        if filename is not None:
            fig.savefig(filename)
            
        if show:
            pyplot.show(fig)
            pyplot.ion()
        
        return fig
        
    def saveKatRun(self, filename):
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
                if "abs:deg" in self.yaxis:
                    out = self.y[:, idx[0]]
                elif "re:im" in self.yaxis:
                    out = self.y[:, idx[0]]
            else: 
                if "abs:deg" in self.yaxis:
                    out = self.y[:, idx[0]] * np.exp(1j*math.pi*self.y[:, idx[1]]/180.0)
                elif "re:im" in self.yaxis :
                    out = self.y[:, idx[0]] + 1j*self.y[:, idx[1]]

            if out is None:
                out = self.y[:, idx]

            if out.size == 1:
                return out[0].squeeze()
            else:
                return out.squeeze()
        else:
            raise  pkex.BasePyKatException("No output by the name '{0}' found in the output".format(str(value)))
            
@canFreeze 
class KatRun2D(object):
    def __init__(self):
        self._unfreeze()
        self.runtime = None
        self.runDateTime = None
        self.startDateTime = datetime.datetime.now()
        self.x = None
        self.y = None
        self.z = None
        self.xlabel = None
        self.ylabel = None
        self.zlabels = None
        self.katScript = None
        self.katVersion = None
        self.stderr = None
        self.stdout = None
        self._freeze()
        
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
    
@canFreeze    
class Signals(object):
    
    @canFreeze 
    class fsig(object):
        def __init__(self, param, name, amplitude, phase, signal):
            self._unfreeze()
            self._params = []
            self.__target = param
            self.__name = name
            self.__amplitude = Param("amp", self, SIfloat(amplitude))
            self.__phase = Param("phase", self, SIfloat(phase))
            self.__removed = False
            self.__signal = signal
            self._freeze()
            
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
                self.__removed = True
        
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
            return self._default_name
        else:
            return self.targets[0].name
            
    @property
    def removed(self): return False # we can never remove the Signal object altogethr just the individual fsig targets

    def remove(self):
        for t in self.targets:
            t.remove()
        
        del self.targets[:]
        
        self.f = None
        
    @property
    def f(self): return self.__f
    @f.setter
    def f(self,value):
        if value is None:
            self.__f.value = None
            return
            
        v = SIfloat(value)
        
        if v <= 0:
            raise pkex.BasePyKatException("Signal frequency must be greater than 0.")
            
        self.__f.value = SIfloat(value)
        
    def __init__(self, kat):
        self._unfreeze()
        self._default_name = "fsignal"
        self.targets = []
        self._params = []
        self.__f = Param("f", self, None)
        self._kat = kat
        self._freeze()
        
    def _register_param(self, param):
        self._params.append(param)
        
    def apply(self, target, amplitude, phase, name=None):
        if target is None:
            raise  pkex.BasePyKatException("No target was specified for signal to be applied")
        
        if name is None:
            name = "sig_" + target._owner().name + "_" + target.name
        
        self.targets.append(Signals.fsig(target, name, amplitude, phase, self))
        
    def getFinesseText(self):
        rtn = []
        
        if self.f.value is not None and self.f is not None:
            if len(self.targets) == 0:
                rtn.append("fsig {name} {frequency}"
                                .format(name = self.name,
                                        frequency=str(self.f.value)))
            else:
                for t in self.targets:
                    rtn.extend(t.getFinesseText())
            
                    rtn.append("fsig {name} {comp} {target} {frequency} {phase} {amplitude}"
                                    .format(name = t.name,
                                            comp=t.owner,
                                            target=t.target,
                                            frequency=str(self.f.value),
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

@canFreeze
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
        self._unfreeze()
        self.__looking = False
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
        self.__variables = {}
        
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
        self.__finesse_dir = None
        
        if kat_code != None and kat_file != None:
            raise pkex.BasePyKatException("Specify either a Kat file or some Kat code, not both.")
        
        if kat_code != None:
            self.parseCommands(kat_code)
        
        if kat_file != None:
            self.loadKatFile(kat_file)
    
        self._freeze()
        
    def deepcopy(self):
        return copy.deepcopy(self)
    
    def getAll(self, type):
        """
        Returns a collection of all objects of the type argument that are
        part of this kat object.
        
        Example:
            # returns all cav commands that are present in this kat object
            cavs = kat.getAll(pykat.commands.cavity)
        """
        items = []
        
        for a in (item for item in self.__class__.__dict__):
            b = getattr(self, a)
            
            if isinstance(b, type):
                items.append(b)

        return tuple(items)
        
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
        values = value.split()
        
        if len(values) == 2:
            scale = values[0]
            mode = values[1]
        else:
            scale = "lin"
            mode = value
        
        if not str(scale) in ["lin", "log"]:
            raise pkex.BasePyKatException("yaxis value '{0}' is not a valid option. Valid options are: lin or log".format(str(mode)))
                
        if not str(mode) in self.yaxis_options:
            raise pkex.BasePyKatException("yaxis value '{0}' is not a valid option. Valid options are: {1}".format(str(mode), ",".join(self.yaxis_options) ))
            
        self.__yaxis = str(scale + " " + mode).strip()

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
    def commands(self):
        return self.__commands.copy()
        
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
        with open(katfile) as f:
            commands= f.read()
            
        self.parseCommands(commands, blocks=blocks)
    
    def parseKatCode(self, code, blocks=None):
        warnings.warn('parseKatCode depreciated, use parseCommands.', stacklevel=2)
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
    
    def getBlocks(self):
        return self.__blocks.keys()
    
    def removeBlock(self, name, failOnBlockNotFound=True):
        
        if name not in self.__blocks:
            if failOnBlockNotFound:
                pkex.PrintError("Error removing block:", pkex.BasePyKatException('Block "{0}" was not found'.format(name)))
                sys.exit(1)
            else:
                return
        
        for o in list(self.__blocks[name].contents):
            self.remove(o)
        
        del self.__blocks[name]
    
    def __str__(self):
         return "".join(self.generateKatScript())
         
    def getVariable(self, name):
        if name not in self.__variables:
            raise pkex.BasePyKatException("Finesse variable `$%s` does not exist." % name)
            
        return self.__variables[name]
    
    def registerVariable(self, name, putter):
        if '$' in name:
            raise pkex.BasePyKatException("Finesse variable name `%s` should not include the `$` symbol as it is added internally." % name)
            
        assert(putter is not None)
        assert(name == putter.name)
        
        if name in self.__variables:
            raise pkex.BasePyKatException("Finesse variable name `%s` already exists." % name)
            
        self.__variables[name] = putter
        
    def unregisterVariable(self, name):
        del self.__variables[name]
    
    def printVariables(self):
        for key in self.__variables:
            print("$" + key, "::::", "owner =", self.__variables[key].owner.name, ", use count =", self.__variables[key].putCount)
    
    def parseCommands(self, commands, blocks=None, addToBlock=None, preserve=False):
        
        commands = str(commands)
        
        try:
            if addToBlock is not None and blocks is not None:
                raise pkex.BasePyKatException("When parsing commands you cannot set both blocks and addToBlock arguments")
            
            # Create a new block if one asked for isn't present
            if addToBlock is not None:
                if addToBlock not in self.__blocks:
                    self.__blocks[addToBlock] = Block(addToBlock)
                
            blockComment = False
        
            commands=self.remove_comments(commands)
        
            commands=self.processConstants(commands)
        
            # Some commands need to be processed after others, and some after that.
            # Here we have two lists of processing priority.
            after_process = ([], [])
        
            for line in commands:
                if len(line.strip()) >= 2:
                    line = line.strip()

                    # Looking for block start or end
                    values = line.split()
                    
                    if addToBlock is None:
                        if values[0] == "%%%":
                            if values[1] == "FTblock":
                                newTag = values[2]
                    
                                if self.__currentTag != None and self.__currentTag != NO_BLOCK: 
                                    warnings.warn("found block {0} before block {1} ended".format(newTag, self.__currentTag))    
                        
                                if newTag in self.__blocks:
                                    #raise pkex.BasePyKatException("Block `{0}` has already been read".format(newTag))
                                    self.__currentTag = newTag
                                    continue # Just add to existing block data
                        
                                self.__blocks[newTag] = Block(newTag) # create new list to store all references to components in block
                                self.__currentTag = newTag
                    
                            if values[1] == "FTend":
                                self.__currentTag = NO_BLOCK
                    
                            continue
                    else:
                        self.__currentTag = addToBlock

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
                    elif(first[:2] == "sq"):
                        obj = pykat.components.squeezer.parseFinesseText(line)
                    elif(first[0:2] == "bs"):
                        obj = pykat.components.beamSplitter.parseFinesseText(line)
                    elif(first[0:2] == "gr"):
                        obj = pykat.components.grating.parseFinesseText(line)
                    elif(first[0:5] == "isol1"):
                        obj = pykat.components.isolator1.parseFinesseText(line)
                    elif(first[0:4] == "isol"):
                        obj = pykat.components.isolator.parseFinesseText(line)
                    elif(first[0:4] == "lens"):
                        obj = pykat.components.lens.parseFinesseText(line)
                    elif(first[0:3] == "mod"):
                        obj = pykat.components.modulator.parseFinesseText(line)
                    elif(first[0:2] == "ad"):
                        obj = pykat.detectors.ad.parseFinesseText(line)
                    elif(first[0:2] == "xd"):
                        obj = pykat.detectors.xd.parseFinesseText(line)
                    elif(first[0:2] == "tf"):
                        obj = pykat.commands.tf.parseFinesseText(line)
                    elif(first[0:2] == "cp"):
                        obj = pykat.detectors.cp.parseFinesseText(line)
                    elif(first[0:2] == "bp"):
                        obj = pykat.detectors.bp.parseFinesseText(line)
                    elif(first[0:4] == "gouy"):
                        obj = pykat.detectors.gouy.parseFinesseText(line)
                    elif(first[0:4] == "beam"):
                        obj = pykat.detectors.beam.parseFinesseText(line)
                    elif(first[0:2] == "pd" and first != "pdtype"):
                        obj = pykat.detectors.pd.parseFinesseText(line)
                    elif(first == "qshot" or first == "qshotS" or first == "qshotN"):
                        obj = pykat.detectors.qshot.parseFinesseText(line)
                    elif(first == "qnoised" or first == "qnoisedS" or first == "qnoisedN"):
                        obj = pykat.detectors.qnoised.parseFinesseText(line)
                    elif(first == "xaxis" or first == "xaxis*"):
                        self.noxaxis = False
                        obj = pykat.commands.xaxis.parseFinesseText(line)
                    elif(first[0:2] == "hd"):
                        obj = pykat.detectors.hd.parseFinesseText(line)
                    elif(first.startswith("qhd")):
                        obj = pykat.detectors.qhd.parseFinesseText(line)
                    elif(first == "x2axis" or first == "x2axis*"):
                        obj = pykat.commands.x2axis.parseFinesseText(line)
                    elif(first == "gauss" or first == "gauss*" or first == "gauss**"):
                        after_process[0].append((line, self.__currentTag))
                    elif(first == "scale"):
                        after_process[1].append((line, self.__currentTag))
                    elif(first == "pdtype"):
                        after_process[0].append((line, self.__currentTag))
                    elif(first == "cav"):
                        after_process[0].append((line, self.__currentTag))
                    elif(first == "func"):
                        after_process[0].append((line, self.__currentTag))
                    elif(first == "variable"):
                        after_process[0].append((line, self.__currentTag))
                    elif(first == "lock"):
                        after_process[0].append((line, self.__currentTag))
                    elif(first == "attr"):
                        after_process[0].append((line, self.__currentTag))
                    elif(first == "noxaxis"):
                        self.noxaxis = True
                    elif(first == "lambda"):
                        v = line.split()
                        self.lambda0 = SIfloat(v[-1])
                    elif(first == "yaxis"):
                        v = line.split(" ", 1)
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
                        after_process[0].append((line, self.__currentTag))
                    elif(first == "noplot"):
                        after_process[1].append((line, self.__currentTag))
                    elif(first == "put" or first == "put*"):
                        after_process[1].append((line, self.__currentTag))
                    else:
                        if self.verbose:
                            print ("Parsing `{0}` into pykat object not implemented yet, added as extra line.".format(line))
                    
                        obj = line
                        # manually add the line to the block contents
                        self.__blocks[self.__currentTag].contents.append(line) 
            
                    if obj != None and not isinstance(obj, six.string_types):
                        if self.hasNamedObject(obj.name):
                            getattr(self, obj.name).remove()
                            if self.verbose:
                                print ("Removed existing object '{0}' of type {1} to add line '{2}'".format(obj.name, obj.__class__, line))

                        self.add(obj, block=self.__currentTag)

            # now process all the varous gauss/attr etc. commands which require
            # components to exist first before they can be processed
            for _ in after_process:
                for item in _:
                    line = item[0]
                    first, rest = line.split(" ",1)
                    block = item[1]
                
                    if first == "gauss" or first == "gauss*" or first == "gauss**":
                        pykat.commands.gauss.parseFinesseText(line, self)
                    
                    elif (first == "cav"):
                        self.add(pykat.commands.cavity.parseFinesseText(line, self), block=block)
                    
                    elif (first == "lock"):
                        self.add(pykat.commands.lock.parseFinesseText(line, self), block=block)
                    
                    elif (first == "func"):
                        self.add(pykat.commands.func.parseFinesseText(line, self), block=block)
                    
                    elif (first == "variable"):
                        self.add(pykat.commands.variable.parseFinesseText(line, self), block=block)
                    
                    elif (first == "noplot"):
                        if not hasattr(self, rest):
                            raise pkex.BasePyKatException("noplot command `{0}` refers to non-existing detector".format(line))
                        
                        getattr(self, rest).noplot = True
                
                    elif (first == "put" or first =="put*"):
                        alt = first == "put*"
                    
                        values = line.split()
                        obj = values[1]
                        target = values[2]
                        variable = values[3]
                    
                        try:
                            if not hasattr(self, obj):
                                raise pkex.BasePyKatException("put command `{0}` refers to non-existing component".format(line))
                    
                            obj = getattr(self, obj)
                    
                            if not hasattr(obj, target):
                                raise pkex.BasePyKatException("put command component `{0}` does not have a parameter `{1}`".format(line, target))
                        
                            target = getattr(obj, target)
                    
                            if not target.isPutable:
                                raise pkex.BasePyKatException("put command `{0}` parameter `{1}` cannot be put to".format(line, target))
                        
                            target.put(self.getVariable(variable.replace('$', '')), alt)
                        
                        except pkex.BasePyKatException as ex:
                            if self.verbose:
                                print("Warning: ", ex.msg)
                                print ("Parsing `{0}` into pykat object not implemented yet, added as extra line.".format(line))
                    
                            obj = line
                            # manually add the line to the block contents
                            self.__blocks[block].contents.append(line)
                        
                 
                    elif (first == "scale"):
                        v = line.split()
                        accepted = ["psd","psd_hf","asd","asd_hf","meter", "ampere", "deg", "rad", "1/deg", "1/rad",]
                
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
                
                        if len(v) == 3:
                            self.signals._default_name = name
                            self.signals.f = SIfloat(v[2])
                        else:
                            if v[2] not in self.__components:
                                raise pkex.BasePyKatException("Could not find the component '{0}'. Line: '{1}'".format(v[2], line))
                
                            comp = self.__components[v[2]]
                
                            if comp._default_fsig() is None:
                                raise pkex.BasePyKatException("Component '{0}' cannot be fsig'd. Line: '{1}'".format(comp.name, line))
                    
                            param_name = None
                            amp = None
                    
                            if len(v) == 3:
                                self.signals._default_name = name
                                freq = SIfloat(v[3])
                            elif len(v) == 5:
                                #param is None
                                freq = SIfloat(v[3])
                                phase = SIfloat(v[4])
                            elif len(v) == 6:
                        
                                try:
                                    SIfloat(v[3])
                                    isFloat = True
                                except:
                                    isFloat = False
                            
                                if isFloat:
                                    freq = SIfloat(v[3])
                                    phase = SIfloat(v[4])
                                    amp = SIfloat(v[5])
                                else:
                                    param_name = v[3]
                                    freq = SIfloat(v[4])
                                    phase = SIfloat(v[5])
                        
                            elif len(v) == 7:
                                param_name = v[3]
                                freq = SIfloat(v[4])
                                phase = SIfloat(v[5])
                                amp = SIfloat(v[6])
                            else:
                                raise pkex.BasePyKatException("'{0}' isnot a valid fsig command".format(line))
                    
                            self.signals.f = freq
                    
                            param = None
                
                            if param_name is None:
                                param = comp._default_fsig()
                            else:
                                for p in comp._params:
                                    if p.canFsig and p.fsigName == param_name:
                                        param = p
                                        break
                    
                                if param is None:
                                    raise pkex.BasePyKatException("Line: '{0}': {1} is not a valid fsig target for {2}".format(line, param_name, comp.name))
                        
                            self.signals.apply(param, amp, phase, name)
                
                    else:
                        raise pkex.BasePyKatException("Haven't handled parsing of '{0}'".format(line))
                    
                self.__currentTag = NO_BLOCK 
        

        except pkex.BasePyKatException as ex:
            pkex.PrintError("Pykat error parsing line: '%s':"%  line, ex)
            sys.exit(1)
            
    def saveScript(self, filename=None):
        """
        Saves the current kat object to a Finesse input file
        """
        try:
            with open(filename,'w') as katfile:
                katScript = "".join(self.generateKatScript())       
                katfile.writelines(katScript)
                katfile.flush()

        except pkex.BasePyKatException as ex:
            print (ex)

    def run(self, plot=None, save_output=False, save_kat=False, kat_name=None, cmd_args=None, getTraceData=False, rethrowExceptions=False, usePipe=True):
        """ 
        Runs the current simulation setup that has been built thus far.
        It returns a KatRun or KatRun2D object which is populated with the various
        data from the simulation run.
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
        
        usePipe           - Version of Finesse 2.1 allow piping between kat and pykat for transfer data an progress.
                            Switching this off means some data and variables won't be populated, but it will work with
                            older versions of Finesse.
        rethrowExceptions - if true exceptions will be thrown again rather than being excepted and calling sys.exit()
        """
        start = time.time()
        
        try:        
            if not hasattr(self, "xaxis") and self.noxaxis != None and self.noxaxis == False:
                raise pkex.BasePyKatException("No xaxis was defined")
            
            if len(self.__katdir) == 0:
                # Get the environment variable for where Finesse is stored
                self.__finesse_dir = os.environ.get('FINESSE_DIR')
                
                if self.__finesse_dir is None :
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
            if self.verbose: print ("Running kat - Started at " + str(datetime.datetime.fromtimestamp(start)))
            
            if hasattr(self, "x2axis") and self.noxaxis == False:
                r = KatRun2D()
            else:
                r = KatRun()
                
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
            if self.__tempname is None:
                katfile = tempfile.NamedTemporaryFile(mode ='w', suffix=".kat", dir=self.__tempdir, delete=False)
            else:
                filepath =os.path.join(self.__tempdir, self.__tempname+".kat" )
                katfile = open( filepath, 'w' ) 
                
            katfile.writelines(r.katScript)
            
            katfile.flush()

            pipe_name = katfile.name + str(uuid.uuid4())
            
            if usePipe:
                cmd=[kat_exec, "--pykat=" + pipe_name]
            else:
                cmd=[kat_exec, "--perl1"]
            
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
            
            if sys.platform == "win32" or sys.platform == "cygwin":
            	# Pipes in windows need to be prefixed with a hidden location.
            	pipe_name = "\\\\.\\pipe\\" + pipe_name

            p = Popen(cmd, stderr=PIPE, stdout=PIPE)

            if self.verbose:
                if self.noxaxis:
                    maxval = 1
                else:
                    maxval = 100
                    
                widgets = [progressbar.Percentage(), ' | ', progressbar.ETA(), ' | ', 'Status']
                
                pb = progressbar.ProgressBar(widgets=widgets, maxval = maxval)

            fifo = None

            _start_kat = time.time()
            
            duration = 2 # Duration for searching for open pipe
            
            try:
                if usePipe == True:
                    while fifo is None:
                        try:
                            if time.time() < _start_kat + duration:
                                time.sleep(0.1)
                                fifo = codecs.open(pipe_name, "r", "utf-8")
                                self.__looking = False
                            else:
                                raise pkex.BasePyKatException("Could not connect to pykat pipe in {0} seconds. Ensure you are using Finesse >= v2.1 and Pykat >= v1.0.0. Or set usePipe=False when making kat object.".format(duration))
                        except FileNotFoundError as ex:
                            if self.verbose:
                                if not self.__looking:
                                    print("Looking for pykat pipe...")
                                    self.__looking = True
                
                if fifo is not None:
                    for line in fifo:
                    
                        #if (sys.version_info < (3, 0)):
                        #    line = line.decode("utf8") # Make sure we're using unicode encoding
                    
                        v = line.split(u":", 1)
                    
                        if len(v) != 2:
                            continue    
                        
                        (tag, line) = v
                    
                        if tag == "version":
                            r.katVersion = line
                        elif tag == "progress" and self.verbose:
                            var = line.split("\t")
                        
                            if len(var) == 3:
                            	pb.currval = int(var[1])
                            	pb.widgets[-1] = var[0] + " " + var[2][:-1]
                            	pb.update()
            finally:
            	if fifo is not None:
            		fifo.close()
			
            (stdout, stderr) = p.communicate()

            r.stdout = stdout.decode('utf-8')
            r.stderr = stderr.decode('utf-8')
            
            k = r.stdout.rfind('computation time:')
            
            if usePipe == False:
                # Set version if not using pipe information
                s = r.stdout.find('(build ') + 7
                e = r.stdout[s:].find(')')
                
                if s == -1 or e == -1:
                    r.katVersion = "Couldn't get version number"
                else:
                    r.katVersion = r.stdout[s:(s+e)]
                
            if k > 0:
                try:
                    line = r.stdout[k:]
                    r.runtime = float(line.split(":")[1].replace("s",""))
                except:
                    r.runtime = 0.0
    
            r.runDateTime = datetime.datetime.now()

            # If Finesse returned an error, just print that and exit!
            if p.returncode != 0:
                raise pkex.FinesseRunError(r.stderr, katfile.name)
            
            self.__prevrunfilename = katfile.name
            
            root = os.path.splitext(katfile.name)
            base = os.path.basename(root[0])
            path = os.path.split(katfile.name)[0]            
            outfile = root[0] + ".out"

            traceData = None
            
            if getTraceData:
                # First see if we have any trace files
                
                traceFiles = [file for file in os.listdir(path) if file.endswith(".trace") and file.startswith(base)]
                
                #print("Found %i trace files" % len(traceFiles))
                #print(path)
                #print(traceFiles)
                
                if len(traceFiles) > 0:
                    import fileinput
                    traceData = []
                    
                    for file in traceFiles:
                        traceData.append({})
                        try:
                            ifile = fileinput.input(os.path.join(path, file))
                    
                            for line in ifile:
                                line = line.strip()
                        
                                if len(line) > 0:
                                    a = line.split(':', 1)
                        
                                    if a[0].isdigit():
                                        values = a[1].split()
                                        
                                        node_name = values[1].split("(")[0]
                                        component_name = values[2].split("(")[0]
                                        
                                        line1x = ifile.readline().replace('(','').replace(')','')
                                        line2x = ifile.readline().replace('(','').replace(')','')
                                        line1y = ifile.readline().replace('(','').replace(')','')
                                        line2y = ifile.readline().replace('(','').replace(')','')
                                        
                                        spqx = line2x.strip().split("gamma")
                                        spqy = line2y.strip().split("gamma")
                                        
                                        qx = spqx[0].split("=")[1].replace('i','j').replace(' ','') 
                                        qy = spqy[0].split("=")[1].replace('i','j').replace(' ','') 
                                        
                                        traceData[-1][node_name] = (pykat.BeamParam(q=complex(qx), wavelength=self.lambda0),
                                                                    pykat.BeamParam(q=complex(qy), wavelength=self.lambda0),
                                                                    component_name)
                                        direc = a[1].split(";")[-1].strip().split(None, 1)[-1]
                                        
                                        traceData[-1][node_name][0].direction = direc
                                        traceData[-1][node_name][1].direction = direc
                            
                        finally:
                            ifile.close()

            
            if save_output:        
                newoutfile = "{0}.out".format(base)
                
                cwd = os.path.os.getcwd()
                newoutfile = os.path.join(cwd,newoutfile)
                
                if os.path.isfile(newoutfile):
                    os.remove(newoutfile)
                    
                os.rename(outfile, newoutfile)

                if self.verbose: print ("\nOutput data saved to '{0}'".format(newoutfile))

            # can't see why this check is needed, causes problems when only detectors
            # not parsed as pykat objects are used
            #if len(self.detectors.keys()) > 0: 
            
            if hasattr(self, "x2axis") and self.noxaxis == False:
                [r.x, r.y, r.z, hdr] = self.readOutFile(outfile)
            
                r.xlabel = hdr[0]
                r.ylabel = hdr[1]
                r.zlabels = [s.strip() for s in hdr[2:]]
                #r.zlabels = map(str.strip, hdr[2:])
            else:
                [r.x, r.y, hdr] = self.readOutFile(outfile)
                
                r.xlabel = hdr[0]
                r.ylabels = [s.strip() for s in hdr[1:]]
                #r.ylabels = map(str.strip, hdr[1:]) // replaced 090415 adf 
                    
            if save_kat:
                if kat_name is None:
                    kat_name = "pykat_output"                
                
                cwd = os.path.os.getcwd()
                newkatfile = os.path.join(cwd, kat_name + ".kat")
                
                if os.path.isfile(newkatfile):
                    os.remove(newkatfile)
                  
                os.rename(katfile.name, newkatfile)         
                
                if self.verbose: print ("Kat file saved to '{0}'".format(newkatfile))
                
            katfile.close()
            perfData = []

            rtn = [r]
            
            if sys.version > '3':
                long = int
            
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
            if rethrowExceptions:
                raise ex 
            else:
                pkex.PrintError("Error from Finesse:", ex)
                
        except pkex.BasePyKatException as ex:
            if rethrowExceptions:
                raise ex 
            else:
                pkex.PrintError("Error from pykat:", ex)
        finally:
            if self.verbose:
                print ("")
                print ("Finished in {0:g} seconds".format(float(time.time() - start)))
            
            
    def remove(self, obj):
        try:
            
            if hasattr(obj, "name") and not isinstance(obj, pykat.finesse.Signals) and not (obj.name in self.__components  or obj.name in self.__detectors or obj.name in self.__commands or obj in self.signals.targets):
                raise pkex.BasePyKatException("'{0}' is not currently in the simulation".format(obj.name))
            
            if hasattr(obj, "removed") and obj.removed:
                raise pkex.BasePyKatException("'{0}' has already been removed".format(obj.name))        

            nodes = None
        
            # store nodes that this componet is attached to as a reference for gui
            if isinstance(obj, Component):
                nodes = self.nodes.getComponentNodes(obj)

            if isinstance(obj, Component):    
                del self.__components[obj.name]
                self.__del_component(obj)
                self.nodes._removeComponent(obj)
                
            elif isinstance(obj, Command):  
                if obj._Command__unique:  
                    del self.__commands[obj.__class__.__name__]
                else:
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
        
            if hasattr(obj, "_on_kat_remove"):
                obj._on_kat_remove()
            
            #import gc
            #print (gc.get_referrers(obj))
            
        except pkex.BasePyKatException as ex:
            pkex.PrintError("Error on removing object:", ex)

    def undumpNodes(self, undumped_name_prefix = "dump"):
        """
        Loops through and removes all dump nodes. Required when running quantum noise
        calculations using qnoised as noise must be injected in where losses occur, such as power
        dumped.
        
        The nodes will be renamed 'dump0', 'dump1', ... If by change a kat file already has a
        node called dump0, dump1, etc. then it will skip that name and move on to the next until
        it finds one that doesn't exist.
        """
        
        i = 0
        node_name = "%s_%i" % (str(undumped_name_prefix), i)
    
        for c in self.components.values():
            for n in c.nodes:
                if n.isDump:
                    while hasattr(self.nodes, node_name):
                        node_name = "%s_%i" % (str(undumped_name_prefix), i)
                        i += 1
                        
                    self.nodes.replaceNode(c, n, self.nodes.createNode(node_name))
        
  
    def getMatrices(self):
        
        import scipy
        from scipy.sparse import coo_matrix
        
        prev = self.noxaxis
        
        self.noxaxis = True
        self.printmatrix = True
        print ("".join(self.generateKatScript()))
        self.verbose = True
        self.run()
        self.printmatrix = None
        self.noxaxis = prev        
        
        if self.__prevrunfilename is None:
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
        
    def add(self, obj, block=NO_BLOCK):
        try:
            self._unfreeze()
            obj.tag = block
            self.__blocks[block].contents.append(obj)
            
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
                
                if obj._Command__unique:
                    self.__commands[obj.__class__.__name__] = obj
                else:
                    self.__commands[obj.name] = obj
                    
                self.__add_command(obj)
                
            else:
                raise pkex.BasePyKatException("Object {0} could not be added".format(obj))
                
            obj._on_kat_add(self)
        
        except pkex.BasePyKatException as ex:
            pkex.PrintError("Error on adding object:", ex)
        finally:
            self._freeze()

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
            
            # ensure we have a shape (num outputs, x, y)
            if len(z.shape) == 2:
                z = z.reshape(1, z.shape[0], z.shape[1])
                
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
            
        if self.yaxis is not None:
            out.append("yaxis {0}\n".format(self.yaxis))
        
        if self.deriv_h is not None:
            out.append("deriv_h {0}\n".format(self.deriv_h))
            
        if self.retrace is not None:
            out.append("retrace {0}\n".format(str(self.retrace)))
            
        if self.printmatrix is not None and self.printmatrix == True:
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
            
            if self.app is None:
                created = True
                self.app = QApplication([""])
                
            if self.pykatgui is None:
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
     
    def _lkat_getProcess(self, callback, **kwargs):
        """
        """
        
        cmd = "\n".join(self.generateKatScript())
        
        return Process(target=f__lkat_process, args=(callback, cmd, kwargs))
    
      
    def _lkat_trace(self, getCavities=True, getNodes=True, getSpaces=True):
        """
        Given the current state of the kat object a new FINESSE process is called and just
        the beam tracing routine is run. The object that is returned contains all the information
        from the beam tracing routine for each node and space components defined as well as cavity 
        commands.   
        """
        if lkat_location is None:
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
            raise pkex.BasePyKatException("Argument is not of type Detector")
        
        name = det.name
        fget = lambda self: self.__get_detector(name)
        
        if hasattr(self, name):
            raise pkex.BasePyKatException("There is something attached to the kat object already called `%s`" % name)
        
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
        
        if com._Command__unique:
            name = com.__class__.__name__
        else:
            name = com.name
            
        if hasattr(self, name):
            raise pkex.BasePyKatException("There is something attached to the kat object already called `%s`" % name)
            
        fget = lambda self: self.__get_command(name)
        
        setattr(self.__class__, name, property(fget))
        setattr(self, '__com_' + name, com)                   

    def __del_command(self, com):

        if not isinstance(com, Command):
            raise exceptions.ValueError("Argument is not of type Command")
        
        if com._Command__unique:
            name = com.__class__.__name__
        else:
            name = com.name
        
        delattr(self.__class__, name)
        delattr(self, '__com_' + name)
        
    def __get_command(self, name):
        return getattr(self, '__com_' + name)            
    
    def __add_component(self, comp):

        if not isinstance(comp, Component):
            raise pkex.BasePyKatException("Argument is not of type Component")
        
        name = comp.name
        
        if hasattr(self, name):
            raise pkex.BasePyKatException("There is something attached to the kat object already called `%s`" % name)
            
        fget = lambda self: self.__get_component(name)
        
        setattr(self.__class__, name, property(fget))
        setattr(self, '__comp_' + name, comp)                   
        
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
