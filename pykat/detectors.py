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
    
    @property
    def name(self): return self.__name        

    def __str__(self): return self.name
    
class photodiode(Detector):

    class __F(list):
        def __init__(self, values=None):
            if values==None:
                values = []
            list.__init__(self,[SIfloat(value) for value in values])
            
    class __Phi(list):
        def __convertValue(self, value):
            if value=="max":
                return value                
            else:
                return SIfloat(value)
                
        def __init__(self, values=None):
            if values==None:
                values = []
            list.__init__(self,[self.__convertValue(value) for value in values])

        def __getitem__(self, key): # probably not needed
            if list.__getitem__(self,key)=="max":
                return list.__getitem__(self,key)
            else:
                return float(list.__getitem__(self,key))
            
    @property
    def f(self): return self.__f

    @property
    def phi(self): return self.__phi
        
    def __init__(self, name, node, senstype="", num_demods=0, demods=[]):        
        Detector.__init__(self, name, node)
        if num_demods>2:
            raise NotImplementedError("pd with more than two demodulations not implemented yet")   
        self.num_demods = num_demods
        self.senstype = senstype

        # every second element into f (starting at 1)
        self.__f = self.__F(demods[::2])
        
        # Every second element into phi (starting at 2)
        self.__phi = self.__Phi(demods[1::2])
        
    @staticmethod
    def parseFinesseText(text): 
        values = text.split(" ")

        if values[0][0:2] != "pd":
            raise exceptions.FinesseParse("'{0}' not a valid photodiode command".format(text))
        if len(values[0])==2:
            __num_demods=0
            __senstype=""
        elif len(values[0])==3 or len(values[0])==4:
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
            __f_phi=range(len(self.f)+len(self.phi))
            __f_phi[::2]=self.f
            __f_phi[1::2]=self.phi
            __f_phi_str = " ".join(map(str, __f_phi))
            
            if self._alternate_beam:
                rtn.append("pd{0}{1} {2} {3} {4}".format(self.senstype, self.num_demods, self.name, __f_phi_str,  self.node.name))
            else:
                rtn.append("pd{0}{1} {2} {3} {4}*".format(self.senstype, self.num_demods, self.name, __f_phi_str,  self.node.name))
            
            if self.noplot:
                rtn.append("noplot {0}".format(self.name))
            
            return rtn
        else:
            return None
    
    def getQGraphicsItem(self):
        if self._svgItem == None:
            self._svgItem = ComponentQGraphicsItem(":/resources/photodiode_red.svg",self,[(-5,11,self.node)])
        
        return self._svgItem    
        
