# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 11:10:01 2013

@author: Daniel
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from pykat import USE_GUI, HAS_OPTIVIS, NoGUIException

import pykat.external.six as six

if six.PY2:
	import exceptions

import warnings
import pykat.exceptions as pkex
import pykat
from pykat.node_network import *
from pykat.exceptions import *
import abc
import copy
from collections import OrderedDict

if HAS_OPTIVIS:
    import optivis.bench.components as optivis_components
    from optivis.view.canvas import OptivisCanvasItemDataType
    from optivis.bench.labels import Label as optivis_label
    from optivis.geometry import Coordinates as optivis_coord
    import PyQt4

from pykat.SIfloat import *
from pykat.param import Param, AttrParam
import weakref
import pykat.exceptions as pkex
from copy import deepcopy
from pykat.freeze import canFreeze

next_component_id = 1

if USE_GUI:
    import pykat.gui.resources
    import pykat.gui.graphics
    from pykat.gui.graphics import *

class NodeGaussSetter(object):
    def __init__(self, component, node):                
        self.__comp = weakref.ref(component)
        self.__node = weakref.ref(node)
        self.__name = None
        
    @property
    def name(self):
        return self.__name
        
    @name.setter
    def name(self, value):
        self.__name = str(value)
    
    @property
    def node(self):
        return self.__node()
    
    @property
    def q(self):
        return self.__node().qx
        
    @q.setter
    def q(self, value):
        self.__node().setGauss(self.__comp(), complex(value))
        
    @property
    def qx(self):
        return self.__node().qx
    @qx.setter
    def qx(self, value):
        self.__node().setGauss(self.__comp(), complex(value))
    
    @property
    def qy(self):
        return self.__node().qy
    @qy.setter
    def qy(self, value):
        self.__node().setGauss(self.__comp(), self.qx, complex(value))
        
id_____pykat_class = 0
 
@canFreeze
class Component(object):
    __metaclass__ = abc.ABCMeta

    def __new__(cls, *args, **kwargs):
        # This creates an instance specific class for the component
        # this enables us to add properties to instances rather than
        # all classes
        global id_____pykat_class
        id_____pykat_class += 1
        cnew_name = str("%s.%s_%i" % (cls.__module__, cls.__name__, id_____pykat_class))
        
        cnew = type(cnew_name, (cls,), {})
        
        o = object.__new__(cnew)
        return o
        
    def __init__(self, name=None):
        self._unfreeze()
        
        self._optivis_component = None
        self.__name = name
        self._svgItem = None
        self._requested_node_names = []
        self._kat = None
        self.tag = None
        self._params = []
        self.__removed = False
        self._default_fsig_param = None
        self.optivisLabelContent = None
        
        # store a unique ID for this component
        global next_component_id
        self.__id = next_component_id
        next_component_id += 1    
         
    def __deepcopy__(self, memo):
        """
        When deep copying a kat object we need to take into account
        the instance specific properties.
        """
        # Here we create a copy of this object based of the base class
        # of this one, otherwise we're making a copy of a copy of a copy...
        result = self.__class__.__new__(self.__class__.__base__)
        result._unfreeze()
        result.__dict__ = copy.deepcopy(self.__dict__, memo)
        
        for _ in result._params:
            _._updateOwner(result)
    
        result._freeze()
        return result
        
    def _register_param(self, param):
        self._params.append(param)
    
    def _default_fsig(self):
        """
        Returns what Finesse internally determines as the default 
        fsig parameter. This is used mainly for parsing fsig command
        lines where no target parameter is stated.
        """
        if self._default_fsig_param != None:
            if not self._default_fsig_param.canFsig:
                raise pkex.BasePyKatException("Default fsig parameter %s is not possible to fsig" % (self.__default_fsig_param.name))
            else:
                return self._default_fsig_param
        else:
            return None
    
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
     
    def _on_kat_remove(self):
        # inform all parameters that we have removed its owner
        # so that it can then warn about any puts/vars/xaxis
        for p in self._params:
            p._onOwnerRemoved()
        
        del self._params[:]

        self.__removed = True
          
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
            name = str(ns.node.name)
            
            if '__nodesetter_' + name in self.__dict__:
                delattr(self, '__nodesetter_' + name)
            
            if name in self.__class__.__dict__:
                delattr(self.__class__, name)
        
        # Now re-add them pointing to the recent nodes
        for node in self.nodes:
            if type(node) != pykat.node_network.DumpNode:
                ns = NodeGaussSetter(self, node)
                self.__add_node_setter(ns)
        
    def __add_node_setter(self, ns):
        self._unfreeze()
        
        if not isinstance(ns, NodeGaussSetter):
            raise exceptions.ValueError("Argument is not of type NodeGaussSetter")
        
        name = str(ns.node.name)
        fget = lambda self: self.__get_node_setter(name)
            
        setattr(self.__class__, name, property(fget))
        setattr(self, '__nodesetter_' + name, ns)                   

        self._freeze()
        
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
        if not USE_GUI:
            raise NoGUIException
            
        return None      
    
    @property
    def removed(self): return self.__removed
    
    @property
    def nodes(self): return self._kat.nodes.getComponentNodes(self) 
    
    @property    
    def name(self): return self.__name      
    
    @property
    def id(self): return self.__id
    
    def __str__(self): return self.name
    
    def remove(self):
        if self.__removed:
            raise pkex.BasePyKatException("{0} has already been marked as removed".format(self.name))
        else:
            self._kat.remove(self)
    
    def getOptivisParameterDict(self):
        if len(self._params) == 0:
            return None
            
        d = OrderedDict()
        
        for p in self._params:
            d[p.name] = OptivisCanvasItemDataType.TEXTBOX
        
        return d
        
    def getOptivisTooltip(self):
        tooltip = "Name: %s" % self.name
        
        for p in self._params:
            if p.value is not None:
                tooltip += "\n%s = %s" %(p.name, str(p.value))
        
        return tooltip

    def setOptivisLabelContent(self):
        """
        Sets default Optivis label contents
        """

        if self.optivisLabelContent is None:
            self.optivisLabelContent = {}

        self.optivisLabelContent["Name"] = self.name


class AbstractMirrorComponent(Component):
    __metaclass__ = abc.ABCMeta
        
    def __init__(self, name, R=None, T=None, L=None, phi=0, Rcx=None, Rcy=None, xbeta=None, ybeta=None, mass=None, r_ap=None, Ix=None, Iy=None, zmech=None, rxmech=None, rymech=None):
        super(AbstractMirrorComponent, self).__init__(name)
 
        if (L != None and R != None and T != None) and SIfloat(R)+SIfloat(T)+SIfloat(L) != 1: 
            raise pkex.BasePyKatException('L+R+T must equal 1 if all are specified at {0}'.format(self.name))
        elif (R != None and L is None and T != None):
            L = 1- (SIfloat(R)+SIfloat(T))
        elif (R is None and L != None and T != None):
            R = 1 - (SIfloat(L)+SIfloat(T))
        elif (R != None and L != None and T is None):
            T = 1 - (SIfloat(L)+SIfloat(R))
        elif (L is None and R is None and T is None):
            raise pkex.BasePyKatException('Must specify at least two of L, R or T')
        
        self.__R = Param("R", self, SIfloat(R))
        self.__T = Param("T", self, SIfloat(T))
        self.__L = Param("L", self, SIfloat(L))
        
        self.__phi = Param("phi", self, SIfloat(phi), canFsig=True, fsig_name="phase")
        self.__Rcx = AttrParam("Rcx", self, SIfloat(Rcx))
        self.__Rcy = AttrParam("Rcy", self, SIfloat(Rcy))
        self.__xbeta = AttrParam("xbeta", self, SIfloat(xbeta), canFsig=True, fsig_name="xbeta")
        self.__ybeta = AttrParam("ybeta", self, SIfloat(ybeta), canFsig=True, fsig_name="ybeta")
        self.__mass = AttrParam("mass", self, SIfloat(mass))
        self.__Ix = AttrParam("Ix", self, SIfloat(Ix))
        self.__Iy = AttrParam("Iy", self, SIfloat(Iy))
        self.__r_ap = AttrParam("r_ap", self, SIfloat(r_ap))
    
        self.__zmech = AttrParam("zmech", self, zmech)
        self.__rxmech = AttrParam("rxmech", self, rxmech)
        self.__rymech = AttrParam("rymech", self, rymech)
        
        self.__z = Param("z", self, None, canFsig=True, isPutable=False, isPutter=False, isTunable=False, fsig_name="z")
        self.__rx = Param("rx", self, None, canFsig=True, isPutable=False, isPutter=False, isTunable=False, fsig_name="rx")
        self.__ry = Param("ry", self, None, canFsig=True, isPutable=False, isPutter=False, isTunable=False, fsig_name="ry")
        self.__Fz = Param("Fz", self, None, canFsig=True, isPutable=False, isPutter=False, isTunable=False, fsig_name="Fz")
        self.__Frx = Param("Frx", self, None, canFsig=True, isPutable=False, isPutter=False, isTunable=False, fsig_name="Frx")
        self.__Fry = Param("Fry", self, None, canFsig=True, isPutable=False, isPutter=False, isTunable=False, fsig_name="Fry")
        
        self.__Fs0 = Param("Fs0", self, None, canFsig=True, isPutable=False, isPutter=False, isTunable=False, fsig_name="Fs0")
        self.__Fs1 = Param("Fs1", self, None, canFsig=True, isPutable=False, isPutter=False, isTunable=False, fsig_name="Fs1")
        self.__Fs0 = Param("s0", self, None, canFsig=True, isPutable=False, isPutter=False, isTunable=False, fsig_name="s0")
        self.__Fs1 = Param("s1", self, None, canFsig=True, isPutable=False, isPutter=False, isTunable=False, fsig_name="s1")
            
        self._default_fsig_param = self.__phi
    
    def setRTL(self, R=None, T=None, L=None):
        if R is not None: self.R.value = R
        if T is not None: self.T.value = T
        if L is not None: self.L.value = L

    def completeRTL(self, R=None, T=None, L=None):
        setValues = sum(x is not None for x in [R,T,L])
        if setValues == 3:
            self.setRTL(R,T,L)
        elif setValues < 2:
            raise pkex.BasePyKatException("must set at least two out of three parameters (R, T, L)")            
        else:
            if R is not None:
                self.R.value = R
            else:
                self.R.value = 1-T-L            
            if T is not None:
                self.T.value = T
            else:
                self.T.value = 1-R-L            
            if L is not None:
                self.L.value = L
            else:
                self.L.value = 1-R-T            

    @property
    def z(self): return self.__z
    @property
    def Fz(self): return self.__Fz
    @property
    def rx(self): return self.__rx
    @property
    def Frx(self): return self.__Frx
    @property
    def ry(self): return self.__ry
    @property
    def Fry(self): return self.__Fry
    
    @property
    def L(self): return self.__L
    @L.setter
    def L(self,value): self.__L.value = SIfloat(value)       
        
    @property
    def r_ap(self): return self.__r_ap
    @r_ap.setter
    def r_ap(self,value): self.__r_ap.value = SIfloat(value)

    @property
    def mass(self): return self.__mass
    @mass.setter
    def mass(self,value): self.__mass.value = SIfloat(value)
    
    @property
    def Ix(self): return self.__Ix
    @Ix.setter
    def Ix(self,value): self.__Ix.value = SIfloat(value)
    
    @property
    def Iy(self): return self.__Iy
    @Iy.setter
    def Iy(self,value): self.__Iy.value = SIfloat(value)
    
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
    def zmech(self): return self.__zmech
    @zmech.setter
    def zmech(self,value): self.__zmech.value = value
    
    @property
    def rxmech(self): return self.__rxmech
    @rxmech.setter
    def rxmech(self,value): self.__rxmech.value = value
    
    @property
    def rymech(self): return self.__rymech
    @rymech.setter
    def rymech(self,value): self.__rymech.value = value
    
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

    def parseAttribute(self, key, value):
        if key in ["Rcx", "Rx", "ROCx", "rx", "rcx", "rocx"]:
            self.Rcx = value
        elif key in ["Rcy", "Ry", "ROCy", "ry", "rcy", "rocy"]:
            self.Rcy = value
        elif key in ["Rc", "ROC", "r", "rc", "roc"]:
            self.Rc = value
        elif key in ["M","m", "Mass", "mass"]:
            self.mass = value
        elif key in ["xbeta", "xBeta", "yaw"]:
            self.xbeta = value
        elif key in ["ybeta", "yBeta", "pitch"]:
            self.ybeta = value
        elif key in ["x_off"]:
            self.x_offset = value
        elif key in ["y_off"]:
            self.y_offset = value
        elif key in ["r_ap"]:
            self.r_ap = value
        elif key in ["Ix","ix"]:
            self.Ix = value
        elif key in ["Iy","iy"]:
            self.Iy = value
        elif key in ["zmech", "mech"]:
            self.zmech = value
        elif key in ["rxmech"]:
            self.rxmech = value
        elif key in ["rymech"]:
            self.rymech = value
        else:
            return False
            
        return True
        
class mirror(AbstractMirrorComponent):
    def __init__(self, name, node1, node2, R=None, T=None, L=None,
                 phi=0, Rcx=None, Rcy=None, xbeta=None, ybeta=None, mass=None, r_ap=None):
        super(mirror, self).__init__(name, R, T, L, phi, Rcx, Rcy, xbeta, ybeta, mass, r_ap)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        self._freeze()
        
    def parseAttributes(self, values):
        
        for key in values.keys():
            if not self.parseAttribute(key, values[key]):
                raise pkex.BasePyKatException("No attribute {0} for mirrors".format(key))
        
    @staticmethod
    def parseFinesseText(text):
        values = text.split()

        if values[0] != "m" and values[0] != "m1" and values[0] != "m2":
            raise pkex.BasePyKatException("'{0}' not a valid Finesse mirror command".format(text))
        
        if len(values) != 7:
            raise pkex.BasePyKatException("Mirror Finesse code format incorrect '{0}'".format(text))

        if len(values[0])==1:
            values.pop(0) # remove initial value
            return mirror(values[0], values[4], values[5], L=None, R=values[1], T=values[2], phi=values[3])
        else:
            if values[0][1]=="1":
                values.pop(0) # remove initial value
                return mirror(values[0], values[4], values[5], R=None, L=values[2], T=values[1], phi=values[3])
            else:
                values.pop(0) # remove initial value
                return mirror(values[0], values[4], values[5], T=None, R=values[1], L=values[2], phi=values[3])

    def getFinesseText(self):
        if abs(self.R + self.T + self.L - 1) > 1e-14:
            raise pkex.BasePyKatException("Mirror {0} has R+T+L = {1}, must equal 1 +- 1e-14".format(self.name, self.R+self.T+self.L))
        
        rtn = []
            
        rtn.append('m {0} {1} {2} {3} {4} {5}'.format(
                self.name, self.R.value, self.T.value, self.phi.value,
                self.nodes[0].name, self.nodes[1].name))

        for p in self._params:
            rtn.extend(p.getFinesseText())
                    
        return rtn
    
    def getOptivisComponent(self):
        self.setOptivisLabelContent()
        
        if self._optivis_component is None:
            self._optivis_component = optivis_components.CavityMirror(name=self.name, aoi=0, tooltip=self.getOptivisTooltip, paramList=self.getOptivisParameterDict(), pykatObject=weakref.ref(self))
            
            lbl = optivis_label(text="", position=optivis_coord(0, -1), item=self._optivis_component)
            lbl.content["Name"] = self.name
            self._optivis_component.labels.append(lbl)
        
        return self._optivis_component
    
    def getOptivisNode(self, mode, kat_node):
        mode = mode.lower()
        
        if mode != "input" and mode.lower() != "output":
            raise pkex.BasePyKatException("Mode must be either input or output not %s" % mode)
        
        if mode == "input":
            if kat_node is self.nodes[0]:
                return self._optivis_component.getInputNode("fr")
            else:
                return self._optivis_component.getInputNode("bk")
                
        elif mode == "output":
            if kat_node is self.nodes[0]:
                return self._optivis_component.getOutputNode("fr")
            else:
                return self._optivis_component.getOutputNode("bk")
          
        
    def getQGraphicsItem(self):
        if not USE_GUI:
            raise NoGUIException
            
        if self._svgItem is None:
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/mirror_flat.svg", self ,[(-4,15,self.nodes[0]), (14,15,self.nodes[1])])
            
        return self._svgItem

class beamSplitter(AbstractMirrorComponent):
    def __init__(self, name, node1, node2, node3, node4, R = None, T = None, L=None, phi = 0, alpha = 0, Rcx = None, Rcy = None, xbeta = None, ybeta = None, mass = None, r_ap = None):
        super(beamSplitter, self).__init__(name, R, T, L, phi, Rcx, Rcy, xbeta, ybeta, mass, r_ap)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        self._requested_node_names.append(node3)
        self._requested_node_names.append(node4)
        self.__alpha = Param("alpha", self, SIfloat(alpha))

        self._freeze()
        
    @property
    def alpha(self): return self.__alpha
    @alpha.setter
    def alpha(self,value): self.__alpha.value = SIfloat(value)
    
    def parseAttributes(self, values):
        
        for key in values.keys():
            if not self.parseAttribute(key, values[key]):
                if key == "alpha":
                    self.alpha = values[key]
                else:
                    raise pkex.BasePyKatException("No attribute {0} for mirrors".format(key))
    
    def getOptivisComponent(self):
        self.setOptivisLabelContent()
        
        if self._optivis_component is None:
            self._optivis_component = optivis_components.BeamSplitter(name=self.name, aoi=-self.alpha, tooltip=self.getOptivisTooltip, paramList=self.getOptivisParameterDict(), pykatObject=weakref.ref(self))
        
        return self._optivis_component
    
    def getOptivisNode(self, mode, kat_node):
        mode = mode.lower()
        
        if mode != "input" and mode.lower() != "output":
            raise pkex.BasePyKatException("Mode must be either input or output")
        
        if mode == "input":
            if kat_node is self.nodes[0]:
                return self._optivis_component.getInputNode("frA")
            elif kat_node is self.nodes[1]:
                return self._optivis_component.getInputNode("frB")
            elif kat_node is self.nodes[2]:
                return self._optivis_component.getInputNode("bkB")
            elif kat_node is self.nodes[3]:
                return self._optivis_component.getInputNode("bkA")
        elif mode == "output":
            if kat_node is self.nodes[0]:
                return self._optivis_component.getOutputNode("frB")
            elif kat_node is self.nodes[1]:
                return self._optivis_component.getOutputNode("frA")
            elif kat_node is self.nodes[2]:
                return self._optivis_component.getOutputNode("bkA")
            elif kat_node is self.nodes[3]:
                return self._optivis_component.getOutputNode("bkB")
                     
    @staticmethod
    def parseFinesseText(text):
        values = text.split()

        if values[0] != "bs" and values[0] != "bs1" and values[0] != "bs2":
            raise pkex.BasePyKatException("'{0}' not a valid Finesse beam splitter command".format(text))
        
        if len(values) != 10:
            raise pkex.BasePyKatException("Beam splitter Finesse code format incorrect '{0}'".format(text))

        if len(values[0])==2:
            values.pop(0) # remove initial value
            return beamSplitter(values[0], values[5], values[6], values[7], values[8],
                                values[1], values[2], None, values[3], values[4])
        elif values[0][2]=="1":
            values.pop(0) # remove initial value
            return beamSplitter(values[0], values[5], values[6], values[7], values[8],
                                None, values[1], values[2], values[3], values[4])
        else:
            values.pop(0) # remove initial value
            return beamSplitter(values[0], values[5], values[6], values[7], values[8],
                                values[1], None, values[2], values[3], values[4])
        
    def getFinesseText(self):
        if abs(self.R + self.T + self.L - 1) > 1e-14:
            raise pkex.BasePyKatException("Beamsplitter {0} has R+T+L = {1}, must equal 1 +- 1e-14".format(self.name, self.R+self.T+self.L))

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
        if not USE_GUI:
            raise NoGUIException
            
        if self._svgItem is None:
            # FIXME: make proper SVG component for beam splitter
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/mirror_flat.svg", self ,[(-4,24,self.nodes[0]), (-4,6,self.nodes[1]), (14,6,self.nodes[2]), (14,24,self.nodes[3])])
            
        return self._svgItem
   
class space(Component):
    def __init__(self, name, node1, node2, L = 0, n = 1, g = None, gx = None, gy = None):
        Component.__init__(self, name)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        self._QItem = None
        self.__L = Param("L", self, SIfloat(L), canFsig=True, fsig_name="phase")
        self.__n = Param("n", self, SIfloat(n))

        self.__gx = AttrParam("gx", self, gx)
        self.__gy = AttrParam("gy", self, gy)
        
        self._default_fsig_param = self.__L

        self._freeze()
        
    @property
    def L(self): return self.__L
    @L.setter
    def L(self,value): self.__L.value = SIfloat(value)
    @property
    def n(self): return self.__n
    @n.setter
    def n(self,value): self.__n.value = SIfloat(value)

    @property
    def g(self):
        if self.__gx.value == self.__gy.value: 
            return self.__g
        else:
            raise pkex.BasePyKatException("Gouy phase in x and y directions are different, use gx and gy properties instead")
            
    @g.setter
    def g(self,value):
        self.__gx.value = SIfloat(value)
        self.__gy.value = SIfloat(value)
        
    @property
    def gx(self): return self.__gx
    @gx.setter
    def gx(self,value): self.__gx.value = SIfloat(value)

    @property
    def gy(self): return self.__gy
    @gy.setter
    def gy(self,value): self.__gy.value = SIfloat(value)
    
    def connectingComponents(self):
        """
        Returns the two components that this space connects.
        """
        a = list(self.nodes[0].components + self.nodes[1].components)
        a = [value for value in a if value != self]
        
        if len(a) != 2:
            raise pkex.BasePyKatException("Space should only connect 2 components")
            
        return a
        
    def parseAttributes(self, values):
        
        for key in values.keys():
            if key in ["gx","gouyx"]:
                self.__gx.value = values[key]
            elif key in ["gy", "gouyy"]:
                self.__gy.value = values[key]
            elif key in ["g","gouy"]:
                self.__gx.value = values[key]
                self.__gy.value = values[key]
            else:
                raise pkex.BasePyKatException("No attribute {0} for spaces".format(key))
                
    @staticmethod
    def parseFinesseText(text):
        values = text.split()

        if values[0] != "s":
            raise pkex.BasePyKatException("'{0}' not a valid Finesse space command".format(text))

        values.pop(0) # remove initial value
        
        if len(values) == 5:
            return space(values[0], values[3], values[4], values[1], values[2])
        elif len(values) == 4:
            return space(values[0], values[2], values[3], values[1])
        else:
            raise pkex.BasePyKatException("Space Finesse code format incorrect '{0}'".format(text))
        
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
        if not USE_GUI:
            raise NoGUIException
            
        if self._QItem is None:
            self._QItem = pykat.gui.graphics.SpaceQGraphicsItem(self)
        
        return self._QItem

class grating(Component):
    def __init__(self, name, node1, node2, node3 = None, node4 = None, n = 2, d = 0, eta_0 = None, eta_1 = None, eta_2 = None, eta_3 = None, rho_0 = None, alpha = None): # TODO: implement Rcx, Rcy and Rc
        Component.__init__(self, name)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)

        if n > 2:
            if node3 != None:
                self._requested_node_names.append(node3)
            else:
                raise pkex.BasePyKatException("Grating node 3 not specified")

        if n > 3:
            if node4 != None:
                self._requested_node_names.append(node4)
            else:
                raise pkex.BasePyKatException("Grating node 4 not specified")

        if n > 4 or n < 2:
            raise pkex.BasePyKatException("Grating must have between 2 and 4 ports")
        
        self.__n = n
        self.__d = Param("d", self, SIfloat(d))
        self.__eta_0 = AttrParam("eta_0", self, SIfloat(eta_0))
        self.__eta_1 = AttrParam("eta_1", self, SIfloat(eta_1))
        self.__eta_2 = AttrParam("eta_2", self, SIfloat(eta_2))
        self.__eta_3 = AttrParam("eta_3", self, SIfloat(eta_3))
        self.__rho_0 = AttrParam("rho_0", self, SIfloat(rho_0))
        self.__alpha = AttrParam("alpha", self, SIfloat(alpha))
        self._svgItem = None

        self._freeze()
        
    @property
    def n(self): return Param('n', self.__n)
    @n.setter
    def n(self, value):
        if value < 2 or value > 4:
            raise pkex.BasePyKatException("Grating must have between 2 and 4 ports")
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
        values = text.split()

        if values[0][0 : 2] != "gr":
            raise pkex.BasePyKatException("'{0}' not a valid Finesse grating command".format(text))

        if len(values[0]) > 2:
            if int(values[0][2]) > 4 or int(values[0][2]) < 2:
                raise pkex.BasePyKatException("Grating must have between 2 and 4 ports")
            else:
                n = int(values[0][2])
        else:
            n = 2

        values.pop(0) # remove initial value
        
        if n == 2:
            if len(values) != 4:
                raise pkex.BasePyKatException("Two port grating must have 2 nodes defined")

            return grating(values[0], values[2], values[3], None, None, n, values[1])
        elif n == 3:
            if len(values) != 5:
                raise pkex.BasePyKatException("Three port grating must have 3 nodes defined")
            
            return grating(values[0], values[2], values[3], values[4], None, n, values[1])
        else:
            if len(values) != 6:
                raise pkex.BasePyKatException("Four port grating must have 4 nodes defined")
            
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
        if not USE_GUI:
            raise NoGUIException
            
        if self._svgItem is None:
            self._svgItem = pykat.gui.graphics.SpaceQGraphicsItem(self) # TODO: make SVG graphic for grating
        
        return self._svgItem

class isolator(Component):
    def __init__(self, name, node1, node2, S = 0, node3="dump", option = 0, L = 0):
        """
        Creates an isolator component. Ligth passes from node 1 to node 2, and from node
        2 to node 3. 

        S         Suppression factor for the reversed direction [power dB].
        L         Loss, fraction of input power loss. Number between 0 and 1.
        option    0: Light passes from node1 to node2, and from node2 to node3. Light
                     is suppressed when going from node3 to node2, and from node2 to
                     node1.
                  1: Light passes from node2 to node1, and from node3 to node2. Light
                     is suppressed when going from node1 to node2, and from node2 to
                     node3. 
        """
        
        Component.__init__(self, name)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        self._requested_node_names.append(node3)
        self._svgItem = None
        self._option = option
        
        self.__S = Param("S",self,SIfloat(S))
        self.__L = Param("L",self,SIfloat(L))

        self._freeze()
        
    @property
    def S(self): return self.__S
    @S.setter
    def S(self, value): self.__S.value = SIfloat(value)

    @property
    def L(self): return self.__L
    @L.setter
    def L(self, value): self.__L.value = SIfloat(value)
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split()

        if values[0] != "isol" and values[0] != "isol*":
            raise pkex.BasePyKatException("'{0}' not a valid Finesse isolator command".format(text))

        if values[0].endswith('*'):
            option = 1
        else:
            option = 0
            
        values.pop(0) # remove initial value
        
        if len(values) == 4:
            return isolator(values[0], values[2], values[3], values[1], option=option)
        elif len(values) == 5:
            # Checking if loss is specified, should be a number.
            if values[2].isnumeric():
                return isolator(values[0], values[3], values[4], values[1],
                                L=values[2], option=option)
            # .. if not a number, it's a node name.
            else:
                return isolator(values[0], values[2], values[3], node3=values[4],
                                S=values[1], option=option)
        elif len(values) == 6:
             return isolator(values[0], values[3], values[4], values[1], values[5], option, values[2])
        else:
            raise pkex.BasePyKatException("Isolator Finesse code format incorrect '{0}'".format(text))
        
    def getFinesseText(self):
        if self._option == 0:
            cmd = "isol"
        elif self._option == 1:
            cmd = "isol*"
            
        rtn = ['{cmd} {0} {1} {2} {3} {4} {5}'.format(self.name, self.S.value, self.L.value, self.nodes[0].name, self.nodes[1].name, self.nodes[2].name, cmd=cmd)]
        
        for p in self._params:
            rtn.extend(p.getFinesseText())

        return rtn
    
    def getOptivisComponent(self):
        self.setOptivisLabelContent()
        
        if self._optivis_component is None:
            self._optivis_component = optivis_components.FaradayIsolator(name=self.name, tooltip=self.getOptivisTooltip, paramList=self.getOptivisParameterDict(), pykatObject=weakref.ref(self))
        
        return self._optivis_component
    
    def getOptivisNode(self, mode, kat_node):
        mode = mode.lower()
        
        if mode != "input" and mode.lower() != "output":
            raise pkex.BasePyKatException("Mode must be either input or output")
        
        if mode == "input":
            if kat_node is self.nodes[0]:
                return self._optivis_component.getInputNode("fr")
            elif kat_node is self.nodes[1]:
                return self._optivis_component.getInputNode("bk")
        elif mode == "output":
            if kat_node is self.nodes[0]:
                return self._optivis_component.getnOutputNode("fr")
            elif kat_node is self.nodes[1]:
                return self._optivis_component.getOutputNode("bk")
     
     
    def getQGraphicsItem(self):
        if not USE_GUI:
            raise NoGUIException
            
        if self._svgItem is None:
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/isolator.svg", self ,[(-4,15,self.nodes[0]), (14,15,self.nodes[1]), (14,24,self.nodes[2])])
        
        return self._svgItem
        
        
        
class isolator1(Component):
    def __init__(self, name, node1, node2, node3, node4):
        """
        Creates a 4-port isolator component.
        """
        
        Component.__init__(self, name)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        self._requested_node_names.append(node3)
        self._requested_node_names.append(node4)
        self._svgItem = None

        self._freeze()
        
    @staticmethod
    def parseFinesseText(text):
        values = text.split()

        if values[0] != "isol1":
            raise pkex.BasePyKatException("'{0}' not a valid Finesse isolator command".format(text))
            
        values.pop(0) # remove initial value
        
        if len(values) == 5:
             return isolator1(values[0], values[1], values[2], values[3], values[4])
        else:
            raise pkex.BasePyKatException("Isolator1 Finesse code format incorrect '{0}'".format(text))
        
    def getFinesseText(self):
        rtn = ['isol1 {0} {1} {2} {3} {4}'.format(self.name, self.nodes[0].name,
                                              self.nodes[1].name, self.nodes[2].name,
                                              self.nodes[3].name)]
        
        for p in self._params:
            rtn.extend(p.getFinesseText())

        return rtn
        
    def getQGraphicsItem(self):
        raise NotImplemented()    
        

class lens(Component):
    def __init__(self, name, node1, node2, f=1, p=None):
        Component.__init__(self, name)
        
        if not ((f is None) ^ (p is None)):
            raise pkex.BasePyKatException("Specify either a focal length or power, not both.")
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        self._svgItem = None
        self.__f = Param("f", self, SIfloat(f))
        self.__p = Param("p", self, SIfloat(p))

        self._freeze()
        
    @property
    def f(self): return self.__f
            
    @f.setter
    def f(self, value):
        self.__f.value = SIfloat(value)
        self.__p.value = None
    
    @property
    def p(self): return self.__p
            
    @p.setter
    def p(self, value):
        self.__p.value = SIfloat(value)
        self.__f.value = None
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split()

        if not values[0].startswith("lens"):
            raise pkex.BasePyKatException("'{0}' not a valid Finesse lens command".format(text))

        alt = values[0].endswith("*")
        
        values.pop(0) # remove initial value
        
        if len(values) == 4:
            if not alt:
                return lens(values[0], values[2], values[3], f=values[1], p=None)
            else:
                return lens(values[0], values[2], values[3], f=None, p=values[1])
        else:
            raise pkex.BasePyKatException("Lens Finesse code format incorrect '{0}'".format(text))
        
    def getFinesseText(self):
        if self.__p.value is None:
            rtn = ['lens {0} {1} {2} {3}'.format(self.name, self.f.value, self.nodes[0].name, self.nodes[1].name)]
        else:
            rtn = ['lens* {0} {1} {2} {3}'.format(self.name, self.p.value, self.nodes[0].name, self.nodes[1].name)]
        
        for p in self._params:
            rtn.extend(p.getFinesseText())
            
        return rtn
    
    def getOptivisComponent(self):
        self.setOptivisLabelContent()
        
        if self._optivis_component is None:
            self._optivis_component = optivis_components.ConvexLens(name=self.name, tooltip=self.getOptivisTooltip, paramList=self.getOptivisParameterDict(), pykatObject=weakref.ref(self))
        
        return self._optivis_component
    
    def getOptivisNode(self, mode, kat_node):
        mode = mode.lower()
        
        if mode != "input" and mode.lower() != "output":
            raise pkex.BasePyKatException("Mode must be either input or output")
        
        if mode == "input":
            if kat_node is self.nodes[0]:
                return self._optivis_component.getInputNode("fr")
            elif kat_node is self.nodes[1]:
                return self._optivis_component.getInputNode("bk")
        elif mode == "output":
            if kat_node is self.nodes[0]:
                return self._optivis_component.getnOutputNode("fr")
            elif kat_node is self.nodes[1]:
                return self._optivis_component.getOutputNode("bk")
                
    def getQGraphicsItem(self):
        if not USE_GUI:
            raise NoGUIException
            
        if self._svgItem is None:
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/lens.svg", self ,[(-4,15,self.nodes[0]), (14,15,self.nodes[1])])
        
        return self._svgItem
        
class modulator(Component):
    def __init__(self, name, node1, node2, f, midx, order, modulation_type='pm', phase=0):
        Component.__init__(self, name)
        
        self._requested_node_names.append(node1)
        self._requested_node_names.append(node2)
        self._svgItem = None
        self.__f = Param("f", self, SIfloat(f), canFsig=True, fsig_name="fre")
        self.__midx = Param("midx", self, SIfloat(midx))
        self.__phase = Param("phase", self, SIfloat(phase), canFsig=True, fsig_name="phase")
        self.__order = 0
        self.order = order
        self.type = modulation_type
        
        self._default_fsig_param = self.__phase
        
        self._freeze()
            
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
        
        try:
            value = int(value)
            
            if value <= 1 and value > 6:
                raise pkex.BasePyKatException("modulator order must be between 1 and 6 or 's' for single sideband")
                
        except ValueError:
            if value != 's' or (isinstance(value, int) and  value <= 1 and value > 6):
                raise pkex.BasePyKatException("modulator order must be between 1 and 6 or 's' for single sideband")

        self.__order = value
        
    
    @property 
    def type(self): return self.__type
    @type.setter
    def type(self, value):
        accepted = ["am", "pm", "yaw", "pitch"]
        
        if value not in accepted:
            raise pkex.BasePyKatException("Modulator type must be: " + ", ".join(accepted))
            
        self.__type = str(value)
    
    @staticmethod
    def parseFinesseText(text):
        v = text.split()

        if v[0] != "mod":
            raise pkex.BasePyKatException("'{0}' not a valid Finesse modulator command".format(text))

        v.pop(0) # remove initial value
        
        if len(v) == 7:
            return modulator(v[0], v[5], v[6], v[1], v[2], v[3], v[4])
        if len(v) == 8:
            return modulator(v[0], v[6], v[7], v[1], v[2], v[3], v[4], phase=v[5])
        else:
            raise pkex.BasePyKatException("Modulator Finesse code format incorrect '{0}'".format(text))
        
    def getFinesseText(self):
        rtn = ['mod {0} {1} {2} {3} {4} {5} {6} {7}'.format(self.name, self.f, self.midx, self.order, self.type, self.phase, self.nodes[0].name, self.nodes[1].name)]
        
        for p in self._params:
            rtn.extend(p.getFinesseText())
            
        return rtn
        
    def getOptivisComponent(self):
        self.setOptivisLabelContent()
        
        if self._optivis_component is None:
            self._optivis_component = optivis_components.ElectroopticModulator(name=self.name, tooltip=self.getOptivisTooltip, paramList=self.getOptivisParameterDict(), pykatObject=weakref.ref(self))
            
        return self._optivis_component
    
    def getOptivisNode(self, mode, kat_node):
        mode = mode.lower()
        
        if mode != "input" and mode.lower() != "output":
            raise pkex.BasePyKatException("Mode must be either input or output")
        
        if mode == "input":
            if kat_node is self.nodes[0]:
                return self._optivis_component.getInputNode("fr")
            elif kat_node is self.nodes[1]:
                return self._optivis_component.getInputNode("bk")
        elif mode == "output":
            if kat_node is self.nodes[0]:
                return self._optivis_component.getnOutputNode("fr")
            elif kat_node is self.nodes[1]:
                return self._optivis_component.getOutputNode("bk")
                
    def getQGraphicsItem(self):
        if not USE_GUI:
            raise NoGUIException
            
        if self._svgItem is None:
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/modulator.svg", self ,[(-4,15,self.nodes[0]), (14,15,self.nodes[1])])
        
        return self._svgItem

class laser(Component):
    def __init__(self, name, node, P=1, f=0, phase=0):
        Component.__init__(self,name)
        
        self._requested_node_names.append(node)
        
        self.__power = Param("P", self, SIfloat(P), canFsig=True, fsig_name="amp")
        self.__f_offset = Param("f", self, f, canFsig=True, fsig_name="freq")
        self.__phase = Param("phase", self, SIfloat(phase), canFsig=True, fsig_name="phase")
        self.__noise = AttrParam("noise", self, None)
        self._svgItem = None
        
        self._default_fsig_param = self.__f_offset

        self._freeze()
        
    @property
    def P(self): return self.__power
    @P.setter
    def P(self,value): self.__power.value = float(value)
    
    @property
    def f(self): return self.__f_offset
    @f.setter
    def f(self,value):
        try:
            self.__f_offset.value = SIfloat(value)
        except:
            self.__f_offset.value = value
    
    @property
    def phase(self): return self.__phase
    @phase.setter
    def phase(self,value): self.__phase.value = float(value)
    
    def parseAttributes(self, values):
        
        for key in values.keys():
            if key == "noise":
                self.__noise.value = values[key]
            else:
                raise pkex.BasePyKatException("No attribute {0} at laser".format(key))
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split()

        if values[0] != "l":
            raise pkex.BasePyKatException("'{0}' not a valid Finesse laser command".format(text))

        values.pop(0) # remove initial value
        
        if len(values) == 5:
            return laser(values[0],values[4],P=values[1],f=values[2],phase=values[3])
        elif len(values) == 4:
            return laser(values[0],values[3],P=values[1],f=values[2], phase=0)
        else:
            raise exceptions.FinesseParse("Laser Finesse code format incorrect '{0}'".format(text))
    
    def getFinesseText(self):
        rtn = ['l {0} {1} {2} {3} {4}'.format(self.name, self.__power.value, self.__f_offset.value, self.__phase.value, self.nodes[0].name)]
        
        for p in self._params:
            rtn.extend(p.getFinesseText())
        
        return rtn

    def getOptivisComponent(self):
        self.setOptivisLabelContent()
        
        if self._optivis_component is None:
            self._optivis_component = optivis_components.Laser(name=self.name, tooltip=self.getOptivisTooltip, paramList=self.getOptivisParameterDict(), pykatObject=weakref.ref(self))
            lbl = optivis_label(text="", position=optivis_coord(0, -1), item=self._optivis_component)
            lbl.content["Name"] = self.name
            self._optivis_component.labels.append(lbl)
            
        return self._optivis_component
    
    def getOptivisNode(self, mode, kat_node):
        mode = mode.lower()
        
        if mode != "input" and mode.lower() != "output":
            raise pkex.BasePyKatException("Mode must be either input or output")
        
        if mode == "input":
            return None
        elif mode == "output":
            return self._optivis_component.getOutputNode("out")
                
    def getQGraphicsItem(self):
        if not USE_GUI:
            raise NoGUIException
            
        if self._svgItem is None:
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/laser.svg", self, [(65,25,self.nodes[0])])
            
        return self._svgItem

class squeezer(Component):
    def __init__(self, name, node, f=0, db=0, angle=0, phase=0, entangled_carrier=False):
        Component.__init__(self,name)
        
        self._requested_node_names.append(node)
        
        self.__f = Param("f", self, 0, canFsig=True, fsig_name="f")
        self.f = f
        
        self.__phase = Param("phase", self, SIfloat(phase), canFsig=True, fsig_name="phase")
        self.__db = Param("db", self, SIfloat(db), canFsig=False, fsig_name="r")
        self.__angle = Param("angle", self, SIfloat(angle), canFsig=False, fsig_name="angle")
        self._svgItem = None
        self.entangled_carrier = entangled_carrier

        self._freeze()
        
        
    @property
    def db(self): return self.__db
    @db.setter
    def db(self,value): self.__db.value = float(value)
    
    @property
    def angle(self): return self.__angle
    @angle.setter
    def angle(self,value): self.__angle.value = float(value)
    
    @property
    def f(self): return self.__f
    @f.setter
    def f(self,value):
        try:
            self.__f.value = SIfloat(value)
        except:
            self.__f.value = value
    
    @property
    def phase(self): return self.__phase
    @phase.setter
    def phase(self,value): self.__phase.value = float(value)
    
    def parseAttributes(self, values):
        pass
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split()
        
        if values[0][:2] != "sq":
            raise pkex.BasePyKatException("'{0}' not a valid Finesse squeezer command".format(text))

        entangled_carrier = values[0].endswith("*")
        
        values.pop(0) # remove initial value
        
        if len(values) == 5:
            return squeezer(values[0], values[4], f=values[1],
                            db=values[2], angle=values[3],
                            entangled_carrier=entangled_carrier)
        else:
            raise exceptions.FinesseParse("Squeezer Finesse code format incorrect '{0}'".format(text))
    
    def getFinesseText(self):
        if self.entangled_carrier:
            rtn = ['sq* {0} {1} {2} {3} {4}'.format(self.name, self.f.value, self.db.value, self.angle.value, self.nodes[0].name)]
        else:
            rtn = ['sq {0} {1} {2} {3} {4}'.format(self.name, self.f.value, self.db.value, self.angle.value, self.nodes[0].name)]
        
        for p in self._params:
            rtn.extend(p.getFinesseText())
        
        return rtn
         
    def getQGraphicsItem(self):
        if not USE_GUI:
            raise NoGUIException
            
        if self._svgItem is None:
            self._svgItem = pykat.gui.graphics.ComponentQGraphicsItem(":/resources/laser.svg", self, [(65,25,self.nodes[0])])
            
        return self._svgItem
            
