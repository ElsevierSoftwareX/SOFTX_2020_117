# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 09:56:53 2013

@author: Daniel
"""
import os
import exceptions
import subprocess
import tempfile
import numpy as np

from colorama import Fore
from pykat.node_network import NodeNetwork
from pykat.detectors import Detector
from pykat.components import Component
from pykat.commands import Command
from pykat.gui.gui import *

class MissingFinesseEnvVar(Exception) :    
    def __str__(self) :
        return "The environment variable FINESSE_DIR was not defined"

class kat:                    
        
    def __init__(self):
    
        self.__components = {}        
        self.__detectors = {}        
        self.__commands = {}        
        self.__gui = None
        self.nodes = NodeNetwork(self)  

        # Various         
        self.phase = None
        self.maxtem = None
        self.noxaxis = None
        
        
    def run(self, printout=1, printerr=1, save_output=False, save_kat=False
            ,kat_name=None) :
        """ Runs the current simulation """
                
        # Get the environment variable for where Finesse is stored
        self.__finesse_dir = os.environ.get('FINESSE_DIR')
                
        if self.__finesse_dir == None :
            raise exceptions.MissingFinesseEnvVar
        
        katfile = tempfile.TemporaryFile(suffix=".kat")
        
        katfile.writelines(self.generate())
        katfile.flush()
        
        kat_exec = os.path.join(self.__finesse_dir,'kat {0}'.format(
                                                            katfile.name))   
                                                            
        p=subprocess.Popen(kat_exec, 
                         stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE)

        [out,err] = p.communicate()
        
        if printout == 1: print Fore.GREEN + out
        if printerr == 1: print Fore.RED + err

        
        [root,ext] = os.path.splitext(katfile.name)
        base = os.path.basename(root)            
        outfile = root + ".out"
            
        [x,y,hdr] = self.readOutFile(outfile)
        
        if save_output:        
            
            newoutfile = "{0}.out".format(base)
            
            cwd = os.path.os.getcwd()
            newoutfile = os.path.join(cwd,newoutfile)
            
            os.rename(outfile, newoutfile)

            print "Output data saved to '{0}'".format(newoutfile)
            
        if save_kat:
            if kat_name == None:
                kat_name = "pykat_output"                
            
            cwd = os.path.os.getcwd()
            newkatfile = os.path.join(cwd, kat_name + ".kat")
            os.rename(katfile.name, newkatfile)         
            
            print "Kat file saved to '{0}'".format(newkatfile)
            

        katfile.close()
        
        return [x,y,hdr]
        
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
            # dont error when adding same command, just replace it
            #if obj.__class__.__name__ in self.__commands :
            #    raise exceptions.ValueError("A command '{0}' has already been added".format([obj.__class__.__name__]))            
            
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
    	y = data[:,1:cols]
        
        return [x, y, hdr]
        
    def generate(self) :
        """ Generates the kat file which can then be run """
        if len(self.__components) == 0 :
            raise exceptions.RuntimeError("No components have been added")

        out = []    
        
        for key in self.__components:       
            txt = self.__components[key].getFinesseText() 
            
            if txt != None:
                if isinstance(txt,list):
                    for t in txt: out.append(t+ "\n")
                else:
                    out.append(txt + "\n")
            
        
        for key in self.__detectors:        
            out.append(self.__detectors[key].getFinesseText() + "\n")
        
        if self.noxaxis != None and self.noxaxis == True:
            out.append("noxaxis\n")

        for key in self.__commands:        
            if self.noxaxis == None or (self.noxaxis == True and isinstance(self.__commands[key], xaxis)):
                out.append(self.__commands[key].getFinesseText() + "\n")
            
        if self.phase != None: out.append("phase {0}\n".format(self.phase))
        if self.maxtem != None: out.append("maxtem {0}\n".format(self.maxtem))            
        
        out.append("gnuterm no\n")
        out.append("pyterm no\n")
        
        return out
        
    def openGUI(self):
        self.__gui = openGUI(self)
    
    def getComponents(self):
        return self.__components.values()
    
    
    def __add_detector(self, det):

        if not isinstance(det, Detector):
            raise exceptions.ValueError("Argument is not of type Command")
        
        name = det.name
        fget = lambda self: self.__get_command(name)
        
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