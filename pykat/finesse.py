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
import subprocess
import tempfile
import numpy as np
import datetime
import pickle
import pykat
import warnings
import re

from pykat.node_network import NodeNetwork
from pykat.detectors import Detector
from pykat.components import Component
from pykat.commands import Command, xaxis
from pykat.gui.gui import pyKatGUI

import pykat.exceptions as pkex

from PyQt4.QtCore import QCoreApplication
from PyQt4.QtGui import QApplication

NO_GUI = False
NO_BLOCK = "NO_BLOCK"

class katRun(object):
    def __init__(self):
        self.runDateTime = datetime.datetime.now()
        self.x = None
        self.y = None
        self.xlabel = None
        self.ylabels = None
        self.katScript = None
        self.katVersion = None
        
    def saveKatRun(self, filename):
        with open(filename,'w') as outfile:
            pickle.dump(self, outfile)
    
    @staticmethod
    def loadKatRun(filename):
        with open(filename,'r') as infile:
            return pickle.load(infile)
    
    def get(self, value): return self[value]
    
    def __getitem__(self, value):
        if value in self.ylabels:
            idx = self.ylabels.index(value)
            if len(self.y.shape) == 1:
                return self.y
            else:
                return self.y[:, idx]
        else:
            raise  pkex.BasePyKatException("No output by the name {0} found".format(value))
      
class katRun2D(object):
    def __init__(self):
        self.runDateTime = datetime.datetime.now()
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
    
    def get(self, value): return self[value]
    
    def __getitem__(self, value):
        idx = [i for i in range(len(self.zlabels)) if self.zlabels[i].split(" ")[0] == str(value)]
        
        if len(idx) > 0:
            return self.z[idx].squeeze()
        else:
            raise  pkex.BasePyKatException("No output by the name {0} found".format(str(value)))
      
class Block:
    def __init__(self, name):
        self.__name = name
        self.contents = [] # List of objects and strings of finesse code
        self.enabled = True 
        
    @property
    def name(self): return self.__name
    
class kat(object):                    
        
    def __init__(self, kat_file=None, kat_code=None, katdir="", katname="", tempdir=None, tempname=None):
        
        self.scene = None # scene object for GUI
        self.verbose = True
        self.__blocks = {} # dictionary of blocks that are used
        self.__components = {}  # dictionary of optical components      
        self.__detectors = {}   # dictionary of detectors
        self.__commands = {}    # dictionary of commands
        self.__extra_lines = [] # an array of strings which are just normal finesse code to include when running
        self.__gui = None
        self.nodes = NodeNetwork(self)  
        self.__katdir = katdir
        self.__katname = katname
        self.__tempdir = tempdir
        self.__tempname = tempname
        self.pykatgui = None
        
        # Various options for running finesse, typicaly the commands with just 1 input
        # and have no name attached to them.
        self.__phase = None
        self.__maxtem = None
        self.__noxaxis = None
        self.__time_code = None
        
        if kat_code != None and kat_file != None:
            raise pkex.BasePyKatException("Specify either a Kat file or some Kat code, not both.")
        
        if kat_code != None:
            self.parseCommands(kat_code)
        
        if kat_file != None:
            self.loadKatFile(kat_file)
        
        cls = type(self)
        self.__class__ = type(cls.__name__, (cls,), {})
        
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

    def logo(self):
        print """                                              ..-
    PyKat                 _                  '(
                          \\`.|\\.__...-\"\"""-_." )
       ..+-----.._        /  ' `            .-'
   . '            `:      7/* _/._\\    \\   (
  (        '::;;+;;:      `-"' =" /,`"" `) /
  L.        \\`:::a:f            c_/     n_'
  ..`--...___`.  .    ,  
   `^-....____:   +."""
    
    def loadKatFile(self, katfile):
        commands=open(katfile).read()
        self.parseCommands(commands)
    
    def parseKatCode(self, code):
        #commands = code.split("\n")
        self.parseCommands(code)
        
    def parseCommands(self, commands):
        blockComment = False

        self.__currentTag= NO_BLOCK
        
        if not (NO_BLOCK in self.__blocks):
            self.__blocks[NO_BLOCK] = Block(NO_BLOCK)
        
        commands=self.remove_comments(commands)
        
        after_process = [] # list of commands that should be processed after 
                           # objects have been set and created
        
        for line in commands.split("\n"):
            #for line in commands:
            if len(line.strip()) >= 2:
                line = line.strip()

                # Looking for block start or end
                values = line.split(" ")
                if values[0] == "%%%":
                    if values[1] == "FTblock":
                        newTag = values[2]
                        
                        if self.__currentTag != None and newTag != self.__currentTag: 
                            warnings.warn("found block {0} before block {1} ended".format(newTag, self.__currentTag))    
                            
                        if newTag in self.__blocks:
                            raise pkex.BasePyKatException("Block `{0}` has already been read")
                            
                        self.__blocks[newTag] = Block(newTag) # create new list to store all references to components in block
                        self.__currentTag = newTag                            
                        
                    if values[1] == "FTend":
                        self.__currentTag = NO_BLOCK
                        
                    continue
                #warnings.warn("current tag {0}".format(self.__currentTag))    

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
                elif(first[0:2] == "pd"):
                    obj = pykat.detectors.photodiode.parseFinesseText(line)
                elif(first == "xaxis" or first == "xaxis*"):
                    obj = pykat.commands.xaxis.parseFinesseText(line)
                elif(first == "x2axis" or first == "x2axis*"):
                    obj = pykat.commands.x2axis.parseFinesseText(line)
                elif(first == "gauss" or first == "gauss*" or first == "gauss**"):
                    after_process.append(line)
                else:
                    if self.verbose:
                        print "Parsing `{0}` into pykat object not implemented yet, added as extra line.".format(line)
                    obj = line
                    # manually add the line to the block contents
                    self.__blocks[self.__currentTag].contents.append(line) 
                
                if obj != None and not isinstance(obj, str):
                    self.add(obj)
                    
        # now process all the varous gauss/attr etc. commands which require
        # components to exist first before they can be processed
        for line in after_process:
            first = line.split(" ",1)[0]
            
            if first == "gauss" or first == "gauss*" or first == "gauss**":
                pykat.commands.gauss.parseFinesseText(line)
            
        self.__currentTag = NO_BLOCK 

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
            print ex
            
    def run(self, printout=0, printerr=0, save_output=False, save_kat=False,kat_name=None) :
        """ 
        Runs the current simulation setup that has been built thus far.
        It returns a katRun or katRun2D object which is populated with the various
        data from the simulation run.
        """
        start = datetime.datetime.now()

        
        try:        
            if not hasattr(self, "xaxis") and self.noxaxis != None and self.noxaxis == False:
                raise pkex.BasePyKatException("No xaxis was defined")
            
            if len(self.__katdir) == 0:
                # Get the environment variable for where Finesse is stored
                self.__finesse_dir = os.environ.get('FINESSE_DIR')
                
                if self.__finesse_dir == None :
                    raise MissingFinesseEnvVar()
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
                raise MissingFinesse()
                
            if self.verbose: print "--------------------------------------------------------------"
            if self.verbose: print "Running kat - Started at " + str(start)
            
            if hasattr(self, "x2axis"):
                r = katRun2D()
            else:
                r = katRun()
                
            r.katScript = "".join(self.generateKatScript())   
            
            # create a kat file which we will write the script into
            if self.__tempname == None:
                katfile = tempfile.NamedTemporaryFile(suffix=".kat", dir=self.__tempdir)
            else:
                filepath =os.path.join(self.__tempdir, self.__tempname+".kat" )
                katfile = open( filepath, 'w' ) 
                
            katfile.writelines(r.katScript)
            katfile.flush()
            
            cmd=[kat_exec, '--perl1']
            
            if self.__time_code:
                cmd.append('--perf-timing')

            cmd.append('--no-backspace')
            # set default format so that less repeated numbers are printed to the
            # output file, should speed up running and parsing of output files
            cmd.append('-format=%.15g')

            cmd.append(katfile.name)
            
            #if self.verbose:
                #print cmd
                
            p=subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            err = ""
            
            if self.verbose: print "Finesse binary output:"
            
            for line in iter(p.stderr.readline, ""):
                #err += line 
                
                if len(line) > 0:
                    if line.rstrip().endswith('%'):
                        vals = line.split("-")
                        action = vals[0].strip()
                        prc = vals[1].strip()[:-1]
                        sys.stdout.write("\r{0} {1}%".format(action, prc))
                    elif line[0] == '*':
                        sys.stdout.write(line)
                        
            [out,errpipe] = p.communicate()
            
            # get the version number
            ix = out.find('build ') + 6
            ix2 = out.find(')',ix)
            r.katVersion = out[ix:ix2]
            
            r.runDateTime = datetime.datetime.now()
            
            if p.returncode != 0:
                raise pkex.FinesseRunError(err, katfile.name)
            
            if printout == 1: print out
            if printerr == 1: print err

            root = os.path.splitext(katfile.name)
            base = os.path.basename(root[0])            
            outfile = root[0] + ".out"
            
            if save_output:        
                newoutfile = "{0}.out".format(base)
                
                cwd = os.path.os.getcwd()
                newoutfile = os.path.join(cwd,newoutfile)
                
                if os.path.isfile(newoutfile):
                    os.remove(newoutfile)
                    
                os.rename(outfile, newoutfile)

                if self.verbose: print "\nOutput data saved to '{0}'".format(newoutfile)
            
            if hasattr(self, "x2axis"):
                [r.x,r.y,r.z,hdr] = self.readOutFile(outfile)
                
                r.xlabel = hdr[0]
                r.ylabel = hdr[1]
                r.zlabels = map(str.strip, hdr[2:])
            else:
                [r.x,r.y,hdr] = self.readOutFile(outfile)
            
                r.xlabel = hdr[0]
                r.ylabels = map(str.strip, hdr[1:])
                            
            if save_kat:
                if kat_name == None:
                    kat_name = "pykat_output"                
                
                cwd = os.path.os.getcwd()
                newkatfile = os.path.join(cwd, kat_name + ".kat")
                
                if os.path.isfile(newkatfile):
                    os.remove(newkatfile)
                  
                os.rename(katfile.name, newkatfile)         
                
                if self.verbose: print "Kat file saved to '{0}'".format(newkatfile)
                

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
            
        except pkex.FinesseRunError as fe:
            print fe
        finally:
            if self.verbose: print ""
            if self.verbose: print "Finished in " + str(datetime.datetime.now()-start)
            
        
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
                
            else :
                raise pkex.BasePyKatException("Object {0} could not be added".format(obj))
                
            obj._on_kat_add(self)
            
        except pkex.BasePyKatException as ex:
            print ex

    def readOutFile(self, filename):
        
        with open(filename,'r') as outfile:
            # read first to lines to get to header line
            outfile.readline()
            outfile.readline()
            
            hdr = outfile.readline().replace('%','').replace('\n','').split(',')
        
        data = np.loadtxt(filename,comments='%',skiprows=4)
        
        if hasattr(self, "x2axis"):
            # need to parse 2D outputs slightly different as they are effectively 2D matrices
            # written in linear form
            x = data[0::(1+self.x2axis.steps),0]
            y = data[0:(1+self.x2axis.steps),1]
            # get rows and columns lined up so that we can reshape a single column of all x/y data
            # into a matrix
            z = data[:,2:].transpose().reshape(data.shape[1]-2, 1+self.xaxis.steps, 1+self.x2axis.steps)
            # once you do this the data for y and x axes need swapping
            z = z.swapaxes(1,2)
            return [x, y, z, hdr]
        else:
            shape_len = len(data.shape)
            
            if shape_len > 1:
                rows,cols = data.shape
                x = data[:,0]
                y = data[:,1:cols].squeeze()
            else:
                rows = 1
                cols = data.shape[0]
                
                x = data[0]
                y = data[1:cols].squeeze()
            
            return [x, y, hdr]
            
    def generateKatScript(self) :
        """ Generates the kat file which can then be run """
        
        out = []    
        
        for key in self.__blocks:
            objs = self.__blocks[key].contents
            
            out.append("%%% FTblock " + key + "\n")
            
            for obj in objs:
                if isinstance(obj, str):
                    out.append(obj + '\n')
                    
                elif isinstance(obj, Component) or isinstance(obj, Detector) or isinstance(obj, Command):
                    txt = obj.getFinesseText() 
                    
                    if txt != None:
                        if isinstance(txt,list):
                            for t in txt:
                                out.append(t + "\n")
                        else:
                            out.append(txt + "\n")
                            
            out.append("%%% FTend " + key + "\n")
        
        if self.noxaxis != None and self.noxaxis == True:
            out.append("noxaxis\n")
            
        # now loop through all the nodes and get any gauss commands
        for key in self.nodes.getNodes():
            txt = self.nodes.getNodes()[key].getFinesseText()
            
            if txt != None:
                if isinstance(txt,list):
                    for t in txt: out.append(t+ "\n")
                else:
                    out.append(txt + "\n")
        
        if self.phase != None: out.append("phase {0}\n".format(self.phase))
        if self.maxtem != None: out.append("maxtem {0}\n".format(self.maxtem))            

        # ensure we don't do any plotting. That should be handled
        # by user themselves
        out.append("gnuterm no\n")
        out.append("pyterm no\n")
        
        return out
        
    def openGUI(self):
        if NO_GUI:
            print  "No PyQt4 module was installed so cannot open a GUI"
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

    def remove_comments(self, string):
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
        return regex.sub(_replacer, string)
