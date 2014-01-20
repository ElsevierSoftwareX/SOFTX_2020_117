# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 11:10:01 2013

@author: Daniel
"""
import exceptions
import pykat.exceptions as pkex
import pykat
from pykat.node_network import *
from pykat.exceptions import *
import abc

import pykat.gui.resources
import pykat.gui.graphics
from pykat.gui.graphics import *
from pykat.SIfloat import *
from pykat.param import Param, AttrParam

next_component_id = 1

class NodeGaussSetter(object):
    def __init__(self, component, node):                
        self.__comp = component
        self.__node = node
    
    @property
    def node(self):
        return self.__node
    
    @property
    def q(self):
        return self.__node.qx
        
    @q.setter
    def q(self, value):
        self.__node.setGauss(self.__comp, complex(value))
        
    @property
    def qx(self):
        return self.__node.qx
    @qx.setter
    def qx(self, value):
        self.__node.setGauss(self.__comp, complex(value))
    
    @property
    def qy(self):
        return self.__node.qy
    @qy.setter
    def qy(self, value):
        self.__node.setGauss(self.__comp, self.qx, complex(value))
        
class Component(object):
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, name):
        self.__name = name
        self._svgItem = None
        self._requested_node_names = []
        self._kat = None
        self.tag = None
        self._params = []
        
        # store a unique ID for this component
        global next_component_id
        self.__id = next_component_id
        next_component_id += 1
        
        # This creates an instance specific class for the component
        # this enables us to add properties to instances rather than
        # all classes
        cls = type(self)
        self.__class__ = type(cls.__name__, (cls,), {})
    
    def _register_param(self, param):
        self._params.append(param)
        
    def _on_kat_add(self, kat):
        """
        Called when this component has been added to a kat object.
        kat is the finesse.kat object which it now belongs to and
        node_array is an array specific to this object which contains
        references to the nodes that are attached to it.
        """
        if self._kat != None:
            raise pkex.BasePyKatException("Component has already been added to a finesse.kat object")
            
        self._kat = kat
        
        kat.nodes.registerComponentNodes(self, self._requested_node_names, self.__on_node_change)
        
    def __on_node_change(self):
        # need to update the node gauss parameter setter members 
        self.__update_node_setters()
        
    def __update_node_setters(self):
        # check if any node setters have already been added. If so we
        # need to remove them. This function should get called if the nodes
        # are updated, either by some function call or the GUI
        key_rm = [k for k in self.__dict__ if k.startswith("__nodesetter_", 0, 13)]
        
        # now we have a list of which to remove
        for key in key_rm:
            ns = self.__dict__[key]
            delattr(self, '__nodesetter_' + ns.node.name)
            delattr(self.__class__, ns.node.name)
        
        for node in self.nodes:
            if type(node) != pykat.node_network.DumpNode:
                ns = NodeGaussSetter(self, node)
                self.__add_node_setter(ns)
        
    def __add_node_setter(self, ns):

        if not isinstance(ns, NodeGaussSetter):
            raise exceptions.ValueError("Argument is not of type NodeGaussSetter")
        
        name = ns.node.name
        fget = lambda self: self.__get_node_setter(name)
        
        setattr(self.__class__, name, property(fget))
        setattr(self, '__nodesetter_' + name, ns)                   

    def __get_node_setter(self, name):
        return getattr(self, '__nodesetter_' + name)   
        
    @staticmethod
    @abc.abstractmethod
    def parseFinesseText(text):
        """Parses Finesse syntax"""
        raise NotImplementedError("This function is not implemented")

    @abc.abstractmethod
    def getFinesseText(self):
        """ Base class for individual Finesse optical components """    
        raise NotImplementedError("This function is not implemented")

    @abc.abstractmethod
    def getQGraphicsItem(self):    
        return None      
    
    @property
    def nodes(self): return self._kat.nodes.getComponentNodes(self) 
    
    @property    
    def name(self): return self.__name      
    
    @property
    def id(self): return self.__id
    
    def __str__(self): return self.name
    
class AbstractMirrorComponent(Component):
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, name, R=0, T=0, phi=0, Rcx=0, Rcy=0, xbeta=0, ybeta=0, mass=0, r_ap=0):
        super(AbstractMirrorComponent, self).__init__(name)

        self.__R = Param("R", self, SIfloat(R))
        self.__T = Param("T", self, SIfloat(T))
        self.__phi = Param("phi", self, SIfloat(phi))
        self.__Rcx = AttrParam("Rcx", self, SIfloat(Rcx))
        self.__Rcy = AttrParam("Rcy", self, SIfloat(Rcy))
        self.__xbeta = AttrParam("xbeta", self, SIfloat(xbeta))
        self.__ybeta = AttrParam("ybeta", self, SIfloat(ybeta))
        self.__mass = AttrParam("mass", self, SIfloat(mass))
        self.__r_ap = AttrParam("r_ap", self, SIfloat(r_ap))
        
    @property
    def r_ap(self): return self.__r_ap
    @r_ap.setter
    def r_ap(self,value): self.__r_ap.value = SIfloat(value)

    @property
    def mass(self): return self.__mass
    @mass.setter
    def mass(self,value): self.__mass.value = SIfloat(value)
    
    @property
    def R(self): return self.__R
    @R.setter
    def R(self,value): self.__R.value = SIfloat(value)
    
    @property
    def T(self): return self.__T
    @T.setter
    def T(self,value): self.__T.value = SIfloat(value)
        
    @property
    def phi(self): return self.__phi
    @phi.setter
    def phi(self,value): self.__phi.value = SIfloat(value)
    
    @property
    def Rcx(self): return self.__Rcx
    @Rcx.setter
    def Rcx(self,value): self.__Rcx.value = SIfloat(value)
    
    @property
    def Rcy(self): return self.__Rcy
    @Rcy.setter
    def Rcy(self,value): self.__Rcy.value = SIfloat(value)
    
    @property
    def xbeta(self): return self.__xbeta
    @xbeta.setter
    def xbeta(self,value): self.__xbeta.value = SIfloat(value)
    
    @property
    def ybeta(self): return self.__ybeta
    @ybeta.setter
    def ybeta(self,value): self.__ybeta.value = SIfloat(value)
    
    @property
    def Rc(self):
        if self.Rcx == self.Rcy:
            return self.Rcx
        else:
            return [self.Rcx, self.Rcy]
    
    @Rc.setter
    def Rc(self,value):
        self.Rcx.value = SIfloat(value)
        self.Rcy.value = SIfloat(value)

class mirror(AbstractMirrorComponent):
    def __init__(self,name,node1,node2,R=0,T=0,phi=0,Rcx=0,Rcy=0,xbeta=0,ybeta=0,mass=0, r_ap=0):
        super(mirror, self).__init__(name, R, T, phi, Rcx, Rcy, xbeta, ybeta, mass, r_ap)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split(" ")

        if values[0] != "m" and values[0] != "m1" and values[0] != "m2":
            raise exceptions.RuntimeError("'{0}' not a valid Finesse mirror command".format(text))
        
        if len(values) != 7:
            raise exceptions.RuntimeError("Mirror Finesse code format incorrect '{0}'".format(text))

        if len(values[0])==1:
            values.pop(0) # remove initial value
            return mirror(values[0], values[4], values[5], R=values[1], T=values[2], phi=values[3])
        else:
            if values[0][1]=="1":
                values.pop(0) # remove initial value
                return mirror(values[0], values[4], values[5], R=1.0-SIfloat(values[1])-SIfloat(values[2]), T=values[1], phi=values[3])
            else:
                values.pop(0) # remove initial value
                return mirror(values[0], values[4], values[5], R=values[1], T=1.0-SIfloat(values[1])-SIfloat(values[2]), phi=values[3])

    def getFinesseText(self):        
        rtn = []
            
        rtn.append('m {0} {1} {2} {3} {4} {5}'.format(
                self.name, self.R.value, self.T.value, self.phi.value,
                self.nodes[0].name, self.nodes[1].name))

        for p in self._params:
            rtn.extend(p.getFinesseText())
                    
        return rtn
        
    def getQGraphicsItem(self):
        if self._svgItem == None:
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/mirror_flat.svg", self ,[(-4,15,self.nodes[0]), (14,15,self.nodes[1])])
            
        return self._svgItem

class beamSplitter(AbstractMirrorComponent):
    def __init__(self, name, node1, node2, node3, node4, R = 0, T = 0, phi = 0, alpha = 0, Rcx = 0, Rcy = 0, xbeta = 0, ybeta = 0, mass = 0, r_ap = 0):
        super(beamSplitter, self).__init__(name, R, T, phi, Rcx, Rcy, xbeta, ybeta, mass, r_ap)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        self._requested_node_names.append(node3)
        self._requested_node_names.append(node4)
             
        self.__alpha = Param("alpha", self, SIfloat(alpha))
        
    @property
    def alpha(self): return self.__alpha
    @alpha.setter
    def alpha(self,value): self.__alpha.value = SIfloat(value)
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split(" ")

        if values[0] != "bs" and values[0] != "bs1" and values[0] != "bs2":
            raise exceptions.RuntimeError("'{0}' not a valid Finesse beam splitter command".format(text))
        
        if len(values) != 10:
            raise exceptions.RuntimeError("Beam splitter Finesse code format incorrect '{0}'".format(text))

        if len(values[0])==2:
            values.pop(0) # remove initial value
            return beamSplitter(values[0], values[5], values[6], values[7], values[8], values[1], values[2], values[3], values[4])
        else:
            if values[0][2]=="1":
                values.pop(0) # remove initial value
                return beamSplitter(values[0], values[5], values[6],
                values[7], values[8], 1.0 - SIfloat(values[1]) -
                SIfloat(values[2]), values[1], values[3], values[4])
            else:
                values.pop(0) # remove initial value
                return beamSplitter(values[0], values[5], values[6],
                values[7], values[8], values[1], 1.0 -
                SIfloat(values[1]) - SIfloat(values[2]), values[3],
                values[4])
            
    def getFinesseText(self):
        rtn = []
            
        rtn.append('bs {0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(
                self.name, self.R.value, self.T.value, self.phi.value,
                self.alpha.value, self.nodes[0].name,
                self.nodes[1].name, self.nodes[2].name,
                self.nodes[3].name))

        for p in self._params:
            rtn.extend(p.getFinesseText())
            
        return rtn
        
    def getQGraphicsItem(self):
        if self._svgItem == None:
            # FIXME: make proper SVG component for beam splitter
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/mirror_flat.svg", self ,[(-4,15,self.nodes[0]), (14,15,self.nodes[1])])
            
        return self._svgItem
   
class space(Component):
    def __init__(self, name, node1, node2, L=0, n=1):
        Component.__init__(self, name)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        self._QItem = None
        self.__L = Param("L", self, SIfloat(L))
        self.__n = Param("n", self, SIfloat(n))
        
    @property
    def L(self): return self.__L
    @L.setter
    def L(self,value): self.__L.value = SIfloat(value)
    @property
    def n(self): return self.__n
    @n.setter
    def n(self,value): self.__n.value = SIfloat(value)
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split(" ")

        if values[0] != "s":
            raise exceptions.RuntimeError("'{0}' not a valid Finesse space command".format(text))

        values.pop(0) # remove initial value
        
        if len(values) == 5:
            return space(values[0], values[3], values[4], values[1], values[2])
        elif len(values) == 4:
            return space(values[0], values[2], values[3], values[1])
        else:
            raise exceptions.RuntimeError("Space Finesse code format incorrect '{0}'".format(text))
        
    def getFinesseText(self):
        rtn = []
        
        if self.__n == 1:
            rtn.append('s {0} {1} {2} {3}'.format(self.name, self.__L.value, self.nodes[0].name, self.nodes[1].name))
        else:
            rtn.append('s {0} {1} {2} {3} {4}'.format(self.name, self.__L.value, self.__n.value, self.nodes[0].name, self.nodes[1].name))
       
        for p in self._params:
            rtn.extend(p.getFinesseText())
            
        return rtn
        
    def getQGraphicsItem(self):
        if self._QItem == None:
            self._QItem = pykat.gui.graphics.SpaceQGraphicsItem(self)
        
        return self._QItem

class grating(Component):
    def __init__(self, name, node1, node2, node3 = None, node4 = None, n = 2, d = 0, eta_0 = 0, eta_1 = 0, eta_2 = 0, eta_3 = 0, rho_0 = 0, alpha = 0): # TODO: implement Rcx, Rcy and Rc
        Component.__init__(self, name)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)

        if n > 2:
            if node3 != None:
                self._requested_node_names.append(node3)
            else:
                raise exceptions.RuntimeError("Grating node 3 not specified")

        if n > 3:
            if node4 != None:
                self._requested_node_names.append(node4)
            else:
                raise exceptions.RuntimeError("Grating node 4 not specified")

        if n > 4 or n < 2:
            raise exceptions.RuntimeError("Grating must have between 2 and 4 ports")
        
        self.__n = n
        self.__d = Param("d", self, SIfloat(d))
        self.__eta_0 = AttrParam("eta_0", self, SIfloat(eta_0))
        self.__eta_1 = AttrParam("eta_1", self, SIfloat(eta_1))
        self.__eta_2 = AttrParam("eta_2", self, SIfloat(eta_2))
        self.__eta_3 = AttrParam("eta_3", self, SIfloat(eta_3))
        self.__rho_0 = AttrParam("rho_0", self, SIfloat(rho_0))
        self.__alpha = AttrParam("alpha", self, SIfloat(alpha))
    
    @property
    def n(self): return Param('n', self.__n)
    @n.setter
    def n(self, value):
        if value < 2 or value > 4:
            raise exceptions.RuntimeError("Grating must have between 2 and 4 ports")
        else:
            self.__n = value
    
    @property
    def d(self): return self.__d
    @d.setter
    def d(self, value): self.__d.value = SIfloat(value)
    
    @property
    def eta_0(self): return self.__eta_0
    @eta_0.setter
    def eta_0(self, value): self.__eta_0.value = SIfloat(value)
    
    @property
    def eta_1(self): return self.__eta_1
    @eta_1.setter
    def eta_1(self, value): self.__eta_1.value = SIfloat(value)
    
    @property
    def eta_2(self): return self.__eta_2
    @eta_2.setter
    def eta_2(self, value): self.__eta_2.value = SIfloat(value)
    
    @property
    def eta_3(self): return self.__eta_3
    @eta_3.setter
    def eta_3(self, value): self.__eta_3.value = SIfloat(value)
    
    @property
    def rho_0(self): return self.__rho_0
    @rho_0.setter
    def rho_0(self, value): self.__rho_0.value = SIfloat(value)
    
    @property
    def alpha(self): return self.__alpha
    @alpha.setter
    def alpha(self, value): self.__alpha.value = SIfloat(value)
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split(" ")

        if values[0][0 : 2] != "gr":
            raise exceptions.RuntimeError("'{0}' not a valid Finesse grating command".format(text))

        if len(values[0]) > 2:
            if int(values[0][2]) > 4 or int(values[0][2]) < 2:
                raise exceptions.RuntimeError("Grating must have between 2 and 4 ports")
            else:
                n = int(values[0][2])
        else:
            n = 2

        values.pop(0) # remove initial value
        
        if n == 2:
            if len(values) != 4:
                raise exceptions.RuntimeError("Two port grating must have 2 nodes defined")

            return grating(values[0], values[2], values[3], None, None, n, values[1])
        elif n == 3:
            if len(values) != 5:
                raise exceptions.RuntimeError("Three port grating must have 3 nodes defined")
            
            return grating(values[0], values[2], values[3], values[4], None, n, values[1])
        else:
            if len(values) != 6:
                raise exceptions.RuntimeError("Four port grating must have 4 nodes defined")
            
            return grating(values[0], values[2], values[3], values[4], values[5], n, values[1])
        
    def getFinesseText(self):
        rtn = []
        
        if self.__n == 2:
            rtn.append('gr{0} {1} {2} {3} {4}'.format(self.__n, self.name, self.__d.value, self.nodes[0].name, self.nodes[1].name))
        elif self.__n == 3:
            rtn.append('gr{0} {1} {2} {3} {4} {5}'.format(self.__n, self.name, self.__d.value, self.nodes[0].name, self.nodes[1].name, self.nodes[2].name))
        else:
            rtn.append('gr{0} {1} {2} {3} {4} {5} {6}'.format(self.__n, self.name, self.__d.value, self.nodes[0].name, self.nodes[1].name, self.nodes[2].name, self.nodes[3].name))
        
        for p in self._params:
            rtn.extend(p.getFinesseText())
        
        return rtn
       
    def getQGraphicsItem(self):
        if self._QItem == None:
            self._QItem = pykat.gui.graphics.SpaceQGraphicsItem(self) # TODO: make SVG graphic for grating
        
        return self._QItem

class isolator(Component):
    def __init__(self, name, node1, node2, S = 0):
        Component.__init__(self, name)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        
        self.__S = Param("S",self,SIfloat(S))
        
    @property
    def S(self): return self.__S
    @S.setter
    def S(self, value): self.__S.value = SIfloat(value)
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split(" ")

        if values[0] != "isol":
            raise exceptions.RuntimeError("'{0}' not a valid Finesse isolator command".format(text))

        values.pop(0) # remove initial value
        
        if len(values) == 4:
            return isolator(values[0], values[2], values[3], values[1])
        else:
            raise exceptions.RuntimeError("Isolator Finesse code format incorrect '{0}'".format(text))
        
    def getFinesseText(self):
        rtn = ['isol {0} {1} {2} {3}'.format(self.name, self.S.value, self.nodes[0].name, self.nodes[1].name)]
        
        for p in self._params:
            rtn.extend(p.getFinesseText())
            
        return rtn

    def getQGraphicsItem(self):
        if self._QItem == None:
            self._QItem = pykat.gui.graphics.SpaceQGraphicsItem(self) # TODO: make isolator graphic
        
        return self._QItem

class lens(Component):
    def __init__(self, name, node1, node2, f):
        Component.__init__(self, name)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        
        self.__f = Param("f", self, SIfloat(f))
        
    @property
    def f(self): return self.__f
    @f.setter
    def f(self, value): self.__f.value = SIfloat(value)
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split(" ")

        if values[0] != "lens":
            raise exceptions.RuntimeError("'{0}' not a valid Finesse lens command".format(text))

        values.pop(0) # remove initial value
        
        if len(values) == 4:
            return lens(values[0], values[2], values[3], values[1])
        else:
            raise exceptions.RuntimeError("Lens Finesse code format incorrect '{0}'".format(text))
        
    def getFinesseText(self):
        rtn = ['lens {0} {1} {2} {3}'.format(self.name, self.f.value, self.nodes[0].name, self.nodes[1].name)]
        
        for p in self._params:
            rtn.extend(p.getFinesseText())
            
        return rtn
        
    def getQGraphicsItem(self):
        if self._QItem == None:
            self._QItem = pykat.gui.graphics.SpaceQGraphicsItem(self) # TODO: make lens graphic
        
        return self._QItem
        
class modulator(Component):
    def __init__(self, name, f, midx, order, modulation_type, node1, node2, phase=0):
        Component.__init__(self, name)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        
        self.__f = Param("f", self, SIfloat(f))
        self.__midx = Param("midx", self, SIfloat(midx))
        self.__phase = Param("phase", self, SIfloat(phase))
        self.__order = order
        self.type = modulation_type
        
            
    @property 
    def f(self): return self.__f
    @f.setter
    def f(self, value): self.__f.value = SIfloat(value)
    
    @property 
    def midx(self): return self.__midx
    @midx.setter
    def midx(self, value): self.__midx.value = SIfloat(value)
    
    @property 
    def phase(self): return self.__phase
    @phase.setter
    def phase(self, value): self.__phase.value = SIfloat(value)
    
    @property 
    def order(self): return self.__order
    @order.setter
    def order(self, value):
        if order <= 1 and order > 6:
            raise pkex.BasePyKatException("modulator order must be between 1 and 6")
            
        self.__order = int(value)
    
    @property 
    def type(self): return self.__type
    @type.setter
    def type(self, value):
        if value != "am" and value != "pm":
            raise pkex.BasePyKatException("Modulator type must be am (amplitude modulation) or pm (phase modulation)")
            
        self.__type = str(value)
    
    @staticmethod
    def parseFinesseText(text):
        v = text.split(" ")

        if v[0] != "mod":
            raise exceptions.RuntimeError("'{0}' not a valid Finesse modulator command".format(text))

        v.pop(0) # remove initial value
        
        if len(v) == 7:
            return modulator(v[0], v[1], v[2], v[3], v[4], v[5], v[6])
        if len(v) == 8:
            return modulator(v[0], v[1], v[2], v[3], v[4], v[6], v[7], phase=v[5])
        else:
            raise pkex.BasePyKatException("Modulator Finesse code format incorrect '{0}'".format(text))
        
    def getFinesseText(self):
        rtn = ['mod {0} {1} {2} {3} {4} {5} {6} {7}'.format(self.name, self.f, self.midx, self.order, self.type, self.phase, self.nodes[0].name, self.nodes[1].name)]
        
        for p in self._params:
            rtn.extend(p.getFinesseText())
            
        return rtn
        
    def getQGraphicsItem(self):
        if self._QItem == None:
            self._QItem = pykat.gui.graphics.SpaceQGraphicsItem(self) # TODO: make lens graphic
        
        return self._QItem

class laser(Component):
    def __init__(self,name,node,P=1,f_offset=0,phase=0):
        Component.__init__(self,name)
        
        self._requested_node_names.append(node)
        
        self.__power = Param("P", self, SIfloat(P))
        self.__f_offset = Param("f", self, SIfloat(f_offset))
        self.__phase = Param("phase", self, SIfloat(phase))
        self.__noise = AttrParam("noise", self, 0)
        
    @property
    def power(self): return self.__power
    @power.setter
    def power(self,value): self.__power.value = float(value)
    
    @property
    def f(self): return self.__f_offset
    @f.setter
    def f(self,value): self.__f_offset.value = float(value)
    
    @property
    def phase(self): return self.__phase
    @phase.setter
    def phase(self,value): self.__phase.value = float(value)
    
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
        rtn = ['l {0} {1} {2} {3} {4}'.format(self.name, self.__power.value, self.__f_offset.value, self.__phase.value, self.nodes[0].name)]
        
        for p in self._params:
            rtn.extend(p.getFinesseText())
        
        return rtn
         
    def getQGraphicsItem(self):
        if self._svgItem == None:
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/laser.svg", self, [(65,25,self.nodes[0])])
            
        return self._svgItem
            
