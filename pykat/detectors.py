# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 09:09:10 2013

@author: Daniel
"""
import exceptions
import pykat.gui.resources

from pykat.utils import *
from pykat.gui.graphics import *
from pykat.node_network import *

class Detector(object) :
    def __init__(self, name,node):
        self.__name = name
        self._svgItem = None
        self._kat = None
        self.noplot = False
        self.enabled = True
        self.tag = None
        self.__node = None
        
        if node.find('*'):
            self._alternate_beam = True
            node.replace('*','')
        
        self.__requested_node = node

    def _on_kat_add(self, kat):
        self.__node = kat.nodes.createNode(self.__requested_node)
    
    @staticmethod
    def parseFinesseText(text):    
        raise NotImplementedError("This function is not implemented")
        
    def getFinesseText(self):
        """ Base class for individual finesse optical components """    
        raise NotImplementedError("This function is not implemented")
        
    def getQGraphicsItem(self):    
        return None
    
    @property 
    def node(self): return self.__node
        
    def __getname(self):
        return self.__name        
        
    name = property(__getname)

class photodiode(Detector):

    class F:
        def __init__(self, values=None):
            if values is None:
                self.values = []
            else:
                self.values = values

        def __len__(self):
            return len(self.values)

        def __getitem__(self, key):
            # if key is of invalid type or value, the list values will raise the error
            return self.values[key]

        def __setitem__(self, key, value):
            print "setting"
            self.values[key] = SIfloat(value)

    
    class Phi():
        def __init__(self, values=None):
            print values
            if values is None:
                self.values = []
            else:
                self.values = values

        def __getitem__(self, key): # probably not needed
            print "boom"
            if self.values[key]=="max":
                return self.values[key]
            else:
                return float(self.values[key])
            
        def __setitem__(self,key,value):
            if value=="max":
                self.values[key] = value
            else:
                self.values[key] = SIfloat(value)
        def append(self, value):
            self.values.append(value)

    @property
    def f(self): return self.F('f', self.__f)
    @f.setter
    def f(self, key, value): self.__f[key]=value

    def __init__(self, name, node, senstype=None, num_demods=0, demods=[]):        
        Detector.__init__(self, name, node)
        if num_demods>2:
            raise NotImplementedError("pd with more than two demodulations not implemented yet")   
        self.num_demods = num_demods
        self.senstype = senstype
        # every second element into f (starting at 1)
        
        self.__f = self.F(demods[::2])
        
        # Every second element into phi (starting at 2)
        self.__phi = self.Phi()
        for i in demods[1::2]:
            self.__phi.append(i)
        print self.__phi
        print self.__phi[0]
        

    @staticmethod
    def parseFinesseText(text): 
        values = text.split(" ")

        if values[0][0:2] != "pd":
            raise exceptions.FinesseParse("'{0}' not a valid photodiode command".format(text))
        if len(values[0])==2:
            __num_demods=0
            __senstype=""
        elif len(values[0])==3 or len(values[0])==4:
            print len(values[0])
            if values[0][2]=="S":
                __senstype="S"
            elif values[0][2]=="N":
                __senstype="N"
            else:
                try:
                    __num_demods=int(values[0][2])
                    __senstype=""
                except ValueError:
                    raise exceptions.FinesseParse("'{0}' not a valid photodiode command".format(text))
            if len(values[0])==4:
                try:
                    __num_demods=int(values[0][3])
                except ValueError:
                    raise exceptions.FinesseParse("'{0}' not a valid photodiode command".format(text))                
        else:
            raise exceptions.FinesseParse("'{0}' not a valid photodiode command".format(text))

        if __num_demods<0 or __num_demods>5:
            raise exceptions.FinesseParse("'{0}' number of demodulations must be >0 and <5".format(text))

        values.pop(0) # remove initial value

        if len(values) == 2 * __num_demods + 1 or len(values) == 2 * __num_demods + 2:
            return photodiode(values[0], values[-1], __senstype, __num_demods, values[1:len(values)-1])
        else:
            raise exceptions.FinesseParse("Photodiode code format incorrect '{0}'".format(text))

        #return photodiode("name", "node", demods)   
        #raise NotImplementedError("This function is not implemented")   
        
    def getFinesseText(self) :
        if self.enabled:    
            rtn = []
            
            if self._alternate_beam:
                rtn.append("pd {0} {1}".format(self.name, self.node.name))
            else:
                rtn.append("pd {0} {1}*".format(self.name, self.node.name))
            
            if self.noplot:
                rtn.append("noplot {0}".format(self.name))
            
            return rtn
        else:
            return None
    
    def getQGraphicsItem(self):
        if self._svgItem == None:
            self._svgItem = ComponentQGraphicsItem(":/resources/photodiode_red.svg",self,[(-5,11,self.node)])
        
        return self._svgItem    
        
