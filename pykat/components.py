# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 11:10:01 2013

@author: Daniel
"""
import exceptions
import pykat
from pykat.node_network import *
from pykat.exceptions import *

import pykat.gui.resources
import pykat.gui.graphics
from pykat.gui.graphics import *
from pykat.SIfloat import *

class Component(object) :
    def __init__(self, name):
        self.__name = name
        self._svgItem = None
        self._nodes = []
        self._requested_node_names = []
        self._kat = None
    
    def _on_kat_add(self, kat):
        """
        Called when this component has been added to a kat object
        """
        self._kat = kat
        
        for node_name in self._requested_node_names:
            self._addNode(node_name)
        
    @staticmethod
    def parseFinesseText(text):    
        raise NotImplementedError("This function is not implemented")
    
    def setAttr(name, value):    
        raise NotImplementedError("This function is not implemented")
        
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
                
            self._nodes.append(n)
        
        return n
        
    def getNodes(self):
        """ Returns a copy of the nodes the component has """
        return self._nodes[:]        
            
    def __getname(self):
        return self.__name      
        
    name = property(__getname)
    
class Param(float):
    def __new__(self,name,value):
        return float.__new__(self,SIfloat(value))
         
    def __init__(self,name,value):
        self.__name = name
        
    name = property(lambda self: self.__name)
           
class mirror(Component):
    def __init__(self,name,node1,node2,R=0,T=0,phi=0,Rcx=0,Rcy=0,xbeta=0,ybeta=0,mass=0, r_ap=0):
        
        Component.__init__(self,name)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)

        self.__r_ap = float(r_ap)        
        self.__mass = float(mass)
        self.__R = float(R)
        self.__T = float(T)
        self.__phi = float(phi)
        self.__Rcx = float(Rcx)
        self.__Rcy = float(Rcy)
        self.__xbeta = float(xbeta)
        self.__ybeta = float(ybeta)
    
    @property
    def r_ap(self): return Param('r_ap', self.__mass)
    @r_ap.setter
    def r_ap(self,value): self.__aperture = float(value)

    @property
    def mass(self): return Param('mass', self.__mass)
    @mass.setter
    def mass(self,value): self.__mass = float(value)
    
    @property
    def R(self): return Param('R', self.__R)
    @R.setter
    def R(self,value): self.__R = float(value)
    
    @property
    def T(self): return Param('T', self.__T)
    @T.setter
    def T(self,value): self.__T = float(value)
        
    @property
    def phi(self): return Param('phi', self.__phi)
    @phi.setter
    def phi(self,value): self.__phi = float(value)
    
    @property
    def Rcx(self): return Param('Rcx', self.__Rcx)
    @Rcx.setter
    def Rcx(self,value): self.__Rcx = float(value)
    
    @property
    def Rcy(self): return Param('Rcy', self.__Rcy)
    @Rcy.setter
    def Rcy(self,value): self.__Rcy = float(value)
    
    @property
    def xbeta(self): return Param('xbeta', self.__xbeta)
    @xbeta.setter
    def xbeta(self,value): self.__xbeta = float(value)
    
    @property
    def ybeta(self): return Param('ybeta', self.__ybeta)
    @ybeta.setter
    def ybeta(self,value): self.__ybeta = float(value)
    
    @property
    def Rc(self):
        if self.Rcx == self.Rcy:
            return self.Rcx
        else:
            return [self.Rcx, self.Rcy]
    
    @Rc.setter
    def Rc(self,value):
        self.Rcx = float(value)
        self.Rcy = float(value)
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split(" ")

        if values[0] != "m":
            raise exceptions.RuntimeError("'{0}' not a valid Finesse mirror command".format(text))

        values.pop(0) # remove initial value
        
        if len(values) != 6:
            raise exceptions.RuntimeError("Mirror Finesse code format incorrect '{0}'".format(text))

        return mirror(values[0], values[4], values[5], R=values[1], T=values[2], phi=values[3])
        
    def getFinesseText(self):        
        rtn = []
        nodes = self.getNodes()
        
        if len(nodes) != 2:
            raise exceptions.RuntimeError("Not enough nodes for mirror")
            
        rtn.append('m {0} {1} {2} {3} {4} {5}'.format(
                self.name, self.__R, self.__T, self.__phi,
                nodes[0].name, nodes[1].name))

        if self.r_ap != 0: rtn.append("attr {0} r_ap {1}".format(self.name,self.__r_ap))            
        if self.mass != 0: rtn.append("attr {0} mass {1}".format(self.name,self.__mass))
        if self.Rcx != 0: rtn.append("attr {0} Rcx {1}".format(self.name,self.__Rcx))
        if self.Rcy != 0: rtn.append("attr {0} Rcy {1}".format(self.name,self.__Rcy))
        if self.xbeta != 0: rtn.append("attr {0} xbeta {1}".format(self.name,self.__xbeta))
        if self.ybeta != 0: rtn.append("attr {0} ybeta {1}".format(self.name,self.__ybeta))
        
        return rtn
        
    def getQGraphicsItem(self):
        if self._svgItem == None:
            nodes = self.getNodes()
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/mirror_flat.svg", self ,[(-4,15,nodes[0]), (14,15,nodes[1])])
            
        return self._svgItem
   
   
class space(Component):
    def __init__(self, name, node1, node2, L=0, n=1):
        Component.__init__(self,name,)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        
        self.__L = SIfloat(L)
        self.__n = SIfloat(n)
        self._QItem = None
        
    @property
    def L(self): return Param('L', self.__L)
    @L.setter
    def L(self,value): self.__L = SIfloat(value)
    @property
    def n(self): return Param('n', self.__n)
    @n.setter
    def n(self,value): self.__n = SIfloat(value)
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split(" ")

        if values[0] != "s":
            raise exceptions.RuntimeError("'{0}' not a valid Finesse space command".format(text))

        values.pop(0) # remove initial value
        
        if len(values) == 5:
            return space(values[0],values[3],values[4],L=values[1],n=values[2])
        elif len(values) == 4:
            return space(values[0],values[2],values[3],L=values[1])
        else:
            raise exceptions.RuntimeError("Space Finesse code format incorrect '{0}'".format(text))
        
    def getFinesseText(self):
        nodes = self.getNodes()
        
        if self.__n == 1:
            return 's {0} {1} {2} {3}'.format(self.name, self.__L, nodes[0].name, nodes[1].name)            
        else:
            return 's {0} {1} {2} {3} {4}'.format(self.name, self.__L, self.__n, nodes[0].name, nodes[1].name)            
       
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
        
        if self._nodes[0] == node_old:
            self._nodes[0] = node_new
        
        if self._nodes[1] == node_old:
            self._nodes[1] = node_new

    
class laser(Component):
    def __init__(self,name,node,P=1,f_offset=0,phase=0):
        Component.__init__(self,name)
        
        self._requested_node_names.append(node)
        
        self.__power = float(P)
        self.__f_offset = float(f_offset)
        self.__phase = float(phase)
        
    @property
    def power(self): return Param('P', self.__power)
    @power.setter
    def power(self,value): self.__power = float(value)
    
    @property
    def f_offset(self): return Param('f', self.__f_offset)
    @f_offset.setter
    def f_offset(self,value): self.__f_offset = float(value)
    
    @property
    def phase(self): return Param('phase', self.__phase)
    @phase.setter
    def phase(self,value): self.__phase = float(value)
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split(" ")

        if values[0] != "l":
            raise exceptions.RuntimeError("'{0}' not a valid Finesse laser command".format(text))

        values.pop(0) # remove initial value
        
        if len(values) == 5:
            return laser(values[0],values[4],P=values[1],f_offset=values[2],phase=values[3])
        elif len(values) == 4:
            return laser(values[0],values[3],P=values[1],f_offset=values[2], phase=0)
        else:
            raise exceptions.FinesseParse("Laser Finesse code format incorrect '{0}'".format(text))
    
    def getFinesseText(self):
        nodes = self.getNodes()
        
        return 'l {0} {1} {2} {3} {4}'.format(self.name, self.__power, self.__f_offset, self.__phase, nodes[0].name)            
         
    def getQGraphicsItem(self):
        if self._svgItem == None:
            nodes = self.getNodes()
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/laser.svg",
                                                   self,[(65,25,nodes[0])])
            
        return self._svgItem
            
