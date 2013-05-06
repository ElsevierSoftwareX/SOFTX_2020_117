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
    def __init__(self, name,node,kat):
        self.__name = name
        self._svgItem = None
        self._kat = kat
        self.noplot = False
        self.enabled = True
        
        kat.add(self)
        
        self.__node = kat.nodes.createNode(node)
        self.__node.connect(self)
        
    def getFinesseText(self):
        """ Base class for individual finesse optical components """    
        raise NotImplementedError("This function is not implemented")
        
    def getQGraphicsItem(self):    
        return None
    
    def getNode(self):
        return self.__node;        
        
    def __getname(self):
        return self.__name        
        
    name = property(__getname)

class photodiode(Detector):
    def __init__(self,kat,name,node) :
        Detector.__init__(self,name,node,kat)
        
        if node.find('*'):
            self._alternate_beam = True
            node.replace('*','')
            
        self.__node = kat.nodes.createNode(node)
                    
                
    def getFinesseText(self) :
        if self.enabled:    
            rtn = []
            
            if self._alternate_beam:
                rtn.append("pd {0} {1}".format(self.name, self.__node.name))
            else:
                rtn.append("pd {0} {1}*".format(self.name, self.__node.name))
            
            if self.noplot:
                rtn.append("noplot {0}".format(self.name))
            
            return rtn
        else:
            return None
    
    def getQGraphicsItem(self):
        if self._svgItem == None:
            self._svgItem = ComponentQGraphicsItem(":/resources/photodiode_red.svg",self,[(-20,0,self.node)])
        
        return self._svgItem    
        