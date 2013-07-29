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
import sys
import os
import exceptions
import subprocess
import tempfile
import numpy as np
import datetime
import pickle

from pykat.node_network import NodeNetwork
from pykat.detectors import Detector
from pykat.components import Component
from pykat.commands import Command, xaxis
from pykat.gui.gui import openGUI

class MissingFinesseEnvVar(Exception) :    
    def __str__(self) :
        return "The environment variable FINESSE_DIR was not defined"

class MissingFinesse(Exception) :    
    def __str__(self) :
        return "Could not find the finesse executable 'kat' in '{0}'," \
               "or you do not have the permissions to run it." \
               .format(os.environ.get('FINESSE_DIR'))
    

class katRun(object):
    def __init__(self):
        self.runDateTime = datetime.datetime.now()
        self.x = None
        self.y = None
        self.xlabel = None
        self.ylabels = None
        self.katScript = None
        self.katVersion = None
        
    def saveKatRun(self, run, filename):
        if not isinstance(run, katRun):
            raise RuntimeError("run object must be a katRun type")
        
        with open(filename,'w') as outfile:
            pickle.dump(run, outfile, pickle.HIGHEST_PROTOCOL)
    
    @staticmethod
    def loadKatRun(filename):
        with open(filename,'r') as infile:
            return pickle.load(infile)
        
        
        
class kat(object):                    
        
    def __init__(self, katexe=""):
        
        self.scene = None # scene object for GUI
        self.__components = {}  # dictionary of optical components      
        self.__detectors = {}   # dictionary of detectors
        self.__commands = {}    # dictionary of commands
        self.__extra_lines = [] # an array of strings which are just normal finesse code to include when running
        self.__gui = None
        self.nodes = NodeNetwork(self)  
        self.__katexe = katexe
        
        # Various         
        self.__phase = None
        self.__maxtem = None
        self.__noxaxis = None
        self.__time_code = None
        
    @property
    def maxtem(self): return self.__maxtem
    @maxtem.setter
    def maxtem(self,value): self.__maxtem = int(value)
    
    @property
    def phase(self): return self.__phase
    @phase.setter
    def phase(self,value): self.__phase = int(value)
    
    @property
    def getPerformanceData(self): return self.__time_code
    @getPerformanceData.setter
    def getPerformanceData(self,value): self.__time_code = bool(value)
    
    @property
    def noxaxis(self): return self.__noxaxis
    @noxaxis.setter
    def noxaxis(self,value): self.__noxaxis = bool(value)
       
    def load(self, katfile):
        """
        Loads the kat file specified which can then be run
        """
        
        with open(katfile) as f:
            for lines in f.readlines():
                self.__extra_lines.append(lines)
        
	   
    def run(self, printout=1, printerr=1, save_output=False, save_kat=False,kat_name=None) :
        """ 
        Runs the current simulation setup that has been built thus far.
        It returns a katRun object which is populated with the various
        data from the simulation run.
        """
        
        r = katRun()
        r.katScript = "".join(self.generateKatScript())       
        
        if len(self.__katexe) == 0:
            # Get the environment variable for where Finesse is stored
            self.__finesse_dir = os.environ.get('FINESSE_DIR')
        
            
            if self.__finesse_dir == None :
                raise MissingFinesseEnvVar()
        
            kat_exec = os.path.join(self.__finesse_dir,'kat.exe') 
            
        else:
            kat_exec = self.__katexe
        
        # check if kat file exists and it is executable by user        
        if not (os.path.isfile(kat_exec) and os.access(kat_exec, os.X_OK)):
            raise MissingFinesse()
        
        # create a kat file which we will write the script into
        katfile = tempfile.TemporaryFile(suffix=".kat")
        katfile.writelines(r.katScript)
        katfile.flush()
        
        flags = "--perl1 "
        
        if self.__time_code:
            flags = flags + " --perf-timing --no-backspace"
            
        kat_exec = "{0} {1} {2}".format(kat_exec, flags, katfile.name)
                                                            
        p=subprocess.Popen(kat_exec, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        for line in iter(p.stderr.readline, ""):
            vals = line.split("-")
            
            if len(vals) == 2:
                action = vals[0].strip()
                prc = vals[1].strip()[:-1]
                
                #sys.stdout.write("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
                sys.stdout.write("\r{0} {1}%".format(action, prc))
                sys.stdout.flush()
            
        [out,err] = p.communicate()
        
        # get the version number
        ix = out.find('build ') + 6
        ix2 = out.find(')',ix)
        r.katVersion = out[ix:ix2]
        
        r.runDateTime = datetime.datetime.now()
        
        # if something has gone wrong, print err regardless
        if p.returncode != 0:
            print err
            return None
        
        if printout == 1: print out
        if printerr == 1: print err

        root = os.path.splitext(katfile.name)
        base = os.path.basename(root[0])            
        outfile = root[0] + ".out"
                
        [r.x,r.y,hdr] = self.readOutFile(outfile)
        
        r.xlabel = hdr[0]
        r.ylabels = hdr[1:]
        
        if save_output:        
            newoutfile = "{0}.out".format(base)
            
            cwd = os.path.os.getcwd()
            newoutfile = os.path.join(cwd,newoutfile)
            
            if os.path.isfile(newoutfile):
                os.remove(newoutfile)
                
            os.rename(outfile, newoutfile)

            print "Output data saved to '{0}'".format(newoutfile)
            
        if save_kat:
            if kat_name == None:
                kat_name = "pykat_output"                
            
            cwd = os.path.os.getcwd()
            newkatfile = os.path.join(cwd, kat_name + ".kat")
            
            if os.path.isfile(newkatfile):
                os.remove(newkatfile)
              
            os.rename(katfile.name, newkatfile)         
            
            print "Kat file saved to '{0}'".format(newkatfile)
            

        katfile.close()
        perfData = []
        
        if self.__time_code:
            perffile = open(root[0] + ".perf",'r')
            
            for l in perffile.readlines():
                vals = l.strip().split(' ')
                perfData.append((vals[0], float(vals[1]), float(vals[2]), float(vals[3])))
                
            return [r, perfData]
        else:
            return r
        
    def add(self, obj) :
        
        if isinstance(obj, Component):
            
            if obj.name in self.__components :
                raise exceptions.ValueError("A component with name '{0}' has already been added".format([obj.name]))            
                        
            self.__components[obj.name] = obj
            self.__add_component(obj)
                        
            
        elif isinstance(obj, Detector):
            
            if obj.name in self.__detectors :
                    raise exceptions.ValueError("A detector '{0}' has already been added".format(obj.name))
                    
            self.__detectors[obj.name] = obj
            self.__add_detector(obj)
            
        elif isinstance(obj, Command):
            
            self.__commands[obj.__class__.__name__] = obj
            self.__add_command(obj)
            
        else :
            raise exceptions.ValueError("Object could not be added")
            
            
            # now we have added the component we need to update the node
            # network
            
       
    
    def readOutFile(self, filename):
        
        outfile = open(filename,'r')
        
        # read first to lines to get to header line
        outfile.readline()
        outfile.readline()
        
        hdr = outfile.readline().replace('%','').replace('\n','').split(',')

        data = np.loadtxt(filename,comments='%')
        rows,cols = data.shape
        
        x = data[:,0]
        y = data[:,1:cols].squeeze()
        
        return [x, y, hdr]
            
    def generateKatScript(self) :
        """ Generates the kat file which can then be run """
        
        out = []    
        
        for key in self.__components:       
            txt = self.__components[key].getFinesseText() 
            
            if txt != None:
                if isinstance(txt,list):
                    for t in txt: out.append(t+ "\n")
                else:
                    out.append(txt + "\n")
            
        
        for key in self.__detectors:
            txt = self.__detectors[key].getFinesseText()
            
            if txt != None:
                if isinstance(txt,list):
                    for t in txt: out.append(t+ "\n")
                else:
                    out.append(txt + "\n")
        
        if self.noxaxis != None and self.noxaxis == True:
            out.append("noxaxis\n")

        for key in self.__commands:        
            if self.noxaxis == None or (self.noxaxis == True and isinstance(self.__commands[key], xaxis)):
                txt = self.__commands[key].getFinesseText()
                
                if txt != None:
                    if isinstance(txt,list):
                        for t in txt: out.append(t+ "\n")
                    else:
                        out.append(txt + "\n")
            
        if self.phase != None: out.append("phase {0}\n".format(self.phase))
        if self.maxtem != None: out.append("maxtem {0}\n".format(self.maxtem))            
        
        # There maybe extra lines we want to include in the kat
        # script which aren't parseable into components, detectors
        # or commands. Typically when something hasn't been fully
        # supported yet. So bung these extra lines on at the end
        for lines in self.__extra_lines:
            out.append(lines)
        
        # ensure we don't do any plotting. That should be handled
        # by user themselves
        out.append("gnuterm no\n")
        out.append("pyterm no\n")
        
        return out
        
    def openGUI(self):
        self.__gui = openGUI(self)
    
    def getComponents(self):
        return self.__components.values()
    
    def hasComponent(self, name):
        return (name in self.__components)
    
    def getNewComponentName(self,prefix):
        '''
        Returns a name for a component which hasn't already been added.
        Returns [prefix] + number, where number is greater than 1. e.g.
        if m1 exists getNewName('m') will return 'm2'
        '''
        n = 1
        name = "{0}{1}".format(prefix, n)
        
        while name in self.__components:
            n += 1
            name = "{0}{1}".format(prefix,n)
        
        return name
    
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
        
    
    def __add_detector(self, det):

        if not isinstance(det, Detector):
            raise exceptions.ValueError("Argument is not of type Detector")
        
        name = det.name
        fget = lambda self: self.__get_detector(name)
        
        setattr(self.__class__, name, property(fget))
        setattr(self, '__det_' + name, det)                   

    def __get_detector(self, name):
        return getattr(self, '__det_' + name) 
        
        
    def __add_command(self, com):

        if not isinstance(com, Command):
            raise exceptions.ValueError("Argument is not of type Command")
        
        name = com.__class__.__name__
        fget = lambda self: self.__get_command(name)
        
        setattr(self.__class__, name, property(fget))
        setattr(self, '__com_' + name, com)                   

    def __get_command(self, name):
        return getattr(self, '__com_' + name)            
    
    def __add_component(self, comp):

        if not isinstance(comp, Component):
            raise exceptions.ValueError("Argument is not of type Component")
            
        fget = lambda self: self.__get_component(comp.name)
        
        setattr(self.__class__, comp.name, property(fget))
        setattr(self, '__comp_' + comp.name, comp)                   

    def __get_component(self, name):
        return getattr(self, '__comp_' + name)        