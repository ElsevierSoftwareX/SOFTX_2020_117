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
        
        if node.find('*'):
            self._alternate_beam = True
            node.replace('*','')
        
        self.__requested_node = node

    def _on_kat_add(self, kat):
        self._node = kat.nodes.createNode(self.__requested_node)
    
    @staticmethod
    def parseFinesseText(text):    
        raise NotImplementedError("This function is not implemented")
        
    def getFinesseText(self):
        """ Base class for individual finesse optical components """    
        raise NotImplementedError("This function is not implemented")
        
    def getQGraphicsItem(self):    
        return None
    
    def getNodes(self):
        return [self._node]
        
    def __getname(self):
        return self.__name        
        
    name = property(__getname)

class photodiode(Detector):
            
    def __init__(self, name, node, type, num_demods, demods):
        
        Detector.__init__(self, name, node)
        if num_demods>2:
            raise NotImplementedError("pd with more than two demodulations not implemented yet")   
        self.num_demods = num_demods
        self.type = type
        self

    @property
    def num_demods(self): return Param('num_demods', self.__num_demods)
    @num_demods.setter
    def num_demods(self,value): self.__num_demods = int(value)
    @property
    def type(self): return Param('type', self.__type)
    @type.setter
    def type(self,value): self.__type = value
    @property
    def f1(self): return Param('f1', self.__f1)
    @f1.setter
    def f1(self,value): self.__f1 = SIfloat(value)
    @property
    def phi1(self): return Param('phi1', self.__phi1)
    @phi1.setter
    def phi1(self,value): self.__phi1 = SIfloat(value)    

    @staticmethod
    def parseFinesseText(text): 
        values = text.split(" ")

        if values[0][0:2] != "pd":
            raise exceptions.FinesseParse("'{0}' not a valid photodiode command".format(text))
        if len(value[0])==2:
            __num_demods=0
            __type=""
        elif len(value[0])==3 or len(value[0])==4:
            if value[0][3]=="S":
                __type="S"
            elif value[0][3]=="N":
                __type="N"
            else:
                try:
                    __num_demods=int(values[0][3])
                    __type=""
                except ValueError:
                    raise exceptions.FinesseParse("'{0}' not a valid photodiode command".format(text))
            if len(value[0])==4:
                try:
                    __num_demods=int(values[0][4])
                except ValueError:
                    raise exceptions.FinesseParse("'{0}' not a valid photodiode command".format(text))                
        else:
            raise exceptions.FinesseParse("'{0}' not a valid photodiode command".format(text))

        if __num_demods<0 or __num_demods>5:
            raise exceptions.FinesseParse("'{0}' number of demodulations must be >0 and <5".format(text))

        values.pop(0) # remove initial value
        
        if len(values) == 2 * __num_demods + 1 or len(values) == 2 * __num_demods + 2:
            return photodiode(value[0], values[-1], __type, __num_demods, values[1:len(values-1)])
        else:
            raise exceptions.FinesseParse("Photodiode code format incorrect '{0}'".format(text))

        #return photodiode("name", "node", demods)   
        #raise NotImplementedError("This function is not implemented")   
        
    def getFinesseText(self) :
        if self.enabled:    
            rtn = []
            
            if self._alternate_beam:
                rtn.append("pd {0} {1}".format(self.name, self._node.name))
            else:
                rtn.append("pd {0} {1}*".format(self.name, self._node.name))
            
            if self.noplot:
                rtn.append("noplot {0}".format(self.name))
            
            return rtn
        else:
            return None
    
    def getQGraphicsItem(self):
        if self._svgItem == None:
            self._svgItem = ComponentQGraphicsItem(":/resources/photodiode_red.svg",self,[(-5,11,self._node)])
        
        return self._svgItem    
        
