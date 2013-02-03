# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 11:10:01 2013

@author: Daniel
"""
import exceptions
import pykat.gui.resources
import pykat

from pykat.gui.graphics import *
from pykat.node_network import *
from PyQt4.QtGui import *
from PyQt4.Qt import *

class Component() :
    def __init__(self, name, kat):
        self.__name = name
        self._svgItem = None
        self.__nodes = []
        self._kat = kat
        
        if not isinstance(kat,pykat.finesse.kat):
            raise exceptions.ValueError("kat argument is not a pykat.finesse.kat object")
            
        kat.add(self)
        
    def getFinesseText(self):
        """ Base class for individual finesse optical components """    
        raise NotImplementedError("This function is not implemented")
        
    def getQGraphicsItem(self):    
        return None      
        
    def _addNode(self, name):
        """ Adds a node in sequential order to the component, i.e. add them
        n1, n2, n3, n4... etc. by the name of the node"""
            
        n = self._kat.nodes.createNode(name)
        
        if n == None:
            raise exceptions.RuntimeError("getNode did not return a node for '{0}'".format(name))
        else:
            n.connect(self)    
            self.__nodes.append(n)
        
        return n
        
    def getNodes(self):
        """ Returns a copy of the nodes the component has """
        return self.__nodes[:]        
            
    def __getname(self):
        return self.__name      
        
    name = property(__getname)
    
class Param:
    def __init__(self,name,value):
        self.value = value
        self.__name = name
    
    def getname(self):
        return self.__name
        
    name = property(getname)
            
   
class mirror(Component):
    def __init__(self,kat,name,node1,node2,R=0,T=0,phi=0,Rcx=0,Rcy=0,xbeta=0,ybeta=0):
        
        Component.__init__(self,name,kat)
        
        self.node1 = self._addNode(node1)
        self.node2 = self._addNode(node2)
        
        self.R = Param('R',R)
        self.T = Param('R',T)
        self.phi = Param('phi',phi)
        self.Rcx = Param('rcx',Rcx)
        self.Rcy = Param('rcy',Rcy)
        self.xbeta = Param('xbeta',xbeta)
        self.ybeta = Param('ybeta',ybeta)
        
    def _getRc(self):
        if self.Rcx == self.Rcy:
            return self.Rcx
        else:
            return [self.Rcx, self.Rcy]
    
    def _setRc(self,value):
        self.Rcx = value
        self.Rcy = value
        
    Rc = property(_getRc,_setRc)
    
    def getFinesseText(self):        
        rtn = []
        rtn.append('m {0} {1} {2} {3} {4} {5}'.format(
                self.name, self.R.value, self.T.value, self.phi.value,
                self.node1.name, self.node2.name))
            
        if self.Rcx != 0: rtn.append("attr {0} Rcx {1}".format(self.name,self.Rcx))
        if self.Rcy != 0: rtn.append("attr {0} Rcy {1}".format(self.name,self.Rcy))
        if self.xbeta != 0: rtn.append("attr {0} xbeta {1}".format(self.name,self.xbeta))
        if self.ybeta != 0: rtn.append("attr {0} ybeta {1}".format(self.name,self.ybeta))
        
        return rtn
        
    def getQGraphicsItem(self):
        if self._svgItem == None:
            self._svgItem = ComponentQGraphicsItem(":/resources/mirror_flat.svg",self
                                                ,[(-20,0,self.node1),(20,0,self.node2)])
        return self._svgItem
   
   
   
class space(Component):
    def __init__(self,kat , name, node1, node2, L=0, n=1):
        Component.__init__(self,name,kat)
        
        self.node1 = self._addNode(node1)
        self.node2 = self._addNode(node2)
        
        self.length = Param('L',L)
        self.refractive_index = Param('n',n)
        
    def getFinesseText(self):
        if self.refractive_index.value == 1:
            return 's {0} {1} {2} {3}'.format(self.name, self.length.value, self.node1.name, self.node2.name)            
        else:
            return 's {0} {1} {2} {3} {4}'.format(self.name, self.length.value, self.refractive_index.value, self.node1.name, self.node2.name)            
       
       
       
       
class laser(Component):
    def __init__(self,kat,name,node,P=1,f_offset=0,phase=0):
        Component.__init__(self,name,kat)
                
        self.node = self._addNode(node)
        
        self.power = Param('P', P)
        self.f_offset = Param('f', f_offset)
        self.phase = Param('phase',phase)
        
    def getFinesseText(self):
        if self.phase.value == 0 :
            return 'l {0} {1} {2} {3}'.format(self.name, self.power.value, self.f_offset.value, self.node.name)            
        else :
            return 'l {0} {1} {2} {4} {3}'.format(self.name, self.power.value, self.f_offset.value, self.phase.value, self.node.name)            
         
    def getQGraphicsItem(self):
        if self._svgItem == None:
            self._svgItem = ComponentQGraphicsItem(":/resources/laser.svg",self,[(70,0,self.node)])
            
        return self._svgItem
            
