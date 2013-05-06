# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 11:10:01 2013

@author: Daniel
"""
import exceptions
import pykat.gui.resources
import pykat
import inspect
import pykat.gui.graphics
from pykat.gui.graphics import *
from pykat.node_network import *
from PyQt4.QtGui import *
from PyQt4.Qt import *

class Component(object) :
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
            raise exceptions.RuntimeError("createNode did not return a node for '{0}'".format(name))
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
    
class Param(float):
    def __new__(self,name,value):
        return float.__new__(self,value)
         
    def __init__(self,name,value):
        self.__name = name
        
    name = property(lambda self: self.__name)
           
class mirror(Component):
    def __init__(self,kat,name,node1,node2,R=0,T=0,phi=0,Rcx=0,Rcy=0,xbeta=0,ybeta=0):
        
        Component.__init__(self,name,kat)
        
        self.node1 = self._addNode(node1)
        self.node2 = self._addNode(node2)
        
        self.__R = R
        self.__T = T
        self.__phi = phi
        self.__Rcx = Rcx
        self.__Rcy = Rcy
        self.__xbeta = xbeta
        self.__ybeta = ybeta
            
    @property
    def R(self):
        return Param('R',self.__R)
    @R.setter
    def R(self,value):
        self.__R = value
    @property
    def T(self): return Param('T', self.__T)
    @T.setter
    def T(self,value): self.__T = value
        
    @property
    def phi(self): return Param('phi', self.__phi)
    @phi.setter
    def phi(self,value): self.__phi = value
    
    @property
    def Rcx(self): return Param('Rcx', self.__Rcx)
    @Rcx.setter
    def Rcx(self,value): self.__Rcx = value
    @property
    def Rcy(self): return Param('Rcy', self.__Rcy)
    @Rcy.setter
    def Rcy(self,value): self.__Rcy = value
    
    
    @property
    def xbeta(self): return Param('xbeta', self.__xbeta)
    @xbeta.setter
    def xbeta(self,value): self.__xbeta = value
    @property
    def ybeta(self): return Param('ybeta', self.__ybeta)
    @ybeta.setter
    def ybeta(self,value): self.__ybeta = value
    
    @property
    def Rc(self):
        if self.Rcx == self.Rcy:
            return self.Rcx
        else:
            return [self.Rcx, self.Rcy]
    
    @Rc.setter
    def Rc(self,value):
        self.Rcx = value
        self.Rcy = value
        
    def getFinesseText(self):        
        rtn = []
        rtn.append('m {0} {1} {2} {3} {4} {5}'.format(
                self.name, self.__R, self.__T, self.__phi,
                self.node1.name, self.node2.name))
            
        if self.Rcx != 0: rtn.append("attr {0} Rcx {1}".format(self.name,self.__Rcx))
        if self.Rcy != 0: rtn.append("attr {0} Rcy {1}".format(self.name,self.__Rcy))
        if self.xbeta != 0: rtn.append("attr {0} xbeta {1}".format(self.name,self.__xbeta))
        if self.ybeta != 0: rtn.append("attr {0} ybeta {1}".format(self.name,self.__ybeta))
        
        return rtn
        
    def getQGraphicsItem(self):
        if self._svgItem == None:
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/mirror_flat.svg",self
                                                ,[(-4,15,self.node1),(14,15,self.node2)])
        return self._svgItem
   
   
class space(Component):
    def __init__(self,kat , name, node1, node2, L=0, n=1):
        Component.__init__(self,name,kat)
        
        self.node1 = self._addNode(node1)
        self.node2 = self._addNode(node2)
        
        self.__L = L
        self.__n = n
        self._QItem = None
        
    @property
    def L(self): return Param('L', self.__L)
    @L.setter
    def L(self,value): self.__L = value
    @property
    def n(self): return Param('n', self.__n)
    @n.setter
    def n(self,value): self.__n = value
    
    def getFinesseText(self):
        if self.__n == 1:
            return 's {0} {1} {2} {3}'.format(self.name, self.__L, self.node1.name, self.node2.name)            
        else:
            return 's {0} {1} {2} {3} {4}'.format(self.name, self.__L, self.__n, self.node1.name, self.node2.name)            
       
    def getQGraphicsItem(self):
        if self._QItem == None:
            self._QItem = pykat.gui.graphics.SpaceQGraphicsItem(self)
        
        return self._QItem  

    def changeNode(self, node_old, node_new):
        '''
        Called when a space's node has been connected
        to another components node
        '''
        node_new.connect(self)
        node_old.disconnect(self)
        
        if self.node1 == node_old:
            self.node1 = node_new
        
        if self.node2 == node_old:
            self.node2 = node_new

    
class laser(Component):
    def __init__(self,kat,name,node,P=1,f_offset=0,phase=0):
        Component.__init__(self,name,kat)
                
        self.node = self._addNode(node)
        
        self.__power = P
        self.__f_offset = f_offset
        self.__phase = phase
        
    @property
    def power(self): return Param('P', self.__power)
    @power.setter
    def power(self,value): self.__power = value
    
    @property
    def f_offset(self): return Param('f', self.__f_offset)
    @f_offset.setter
    def f_offset(self,value): self.__f_offset = value
    
    @property
    def phase(self): return Param('phase', self.__phase)
    @phase.setter
    def phase(self,value): self.__phase = value
    
    def getFinesseText(self):
        if self.__phase == 0 :
            return 'l {0} {1} {2} {3}'.format(self.name, self.__power, self.__f_offset, self.node.name)            
        else :
            return 'l {0} {1} {2} {4} {3}'.format(self.name, self.__power, self.__f_offset, self.__phase, self.node.name)            
         
    def getQGraphicsItem(self):
        if self._svgItem == None:
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/laser.svg",
                                                   self,[(65,25,self.node)])
            
        return self._svgItem
            
