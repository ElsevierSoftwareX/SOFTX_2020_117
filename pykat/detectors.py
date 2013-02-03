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
from PyQt4.QtGui import *
from PyQt4.Qt import *

class Detector() :
    def __init__(self, name,node,kat):
        self.__name = name
        self._svgItem = None
        self._kat = kat
        
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
        if self._alternate_beam:
            return "pd {0} {1}".format(self.name, self.__node.name)
        else:
            return "pd {0} {1}*".format(self.name, self.__node.name)
            
    def getQGraphicsItem(self):
        if self._svgItem == None:
            self._svgItem = ComponentQGraphicsItem(":/resources/photodiode_red.svg",self,[(-20,0,self.node)])
        
        return self._svgItem    
        