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
from pykat.param import Param


class Detector(object) :
    def __init__(self, name,node):
        self.__name = name
        self._svgItem = None
        self._kat = None
        self.noplot = False
        self.enabled = True
        self.tag = None
        self.__node = None
        self._params = []
        self._mask = {}
        self.__scale = ""
        
        if node.find('*'):
            self._alternate_beam = True
            node.replace('*','')
        
        self.__requested_node = node
    
    def _register_param(self, param):
        self._params.append(param)
        
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
    def scale(self): return self.__scale
    @scale.setter
    def scale(self, value):
        self.__scale = value

    @property 
    def node(self): return self.__node
    @node.setter
    def node(self, value):
        if value in self._kat.nodes:
            self.__node = self._kat.nodes[value]
        else:
            raise pkex.BasePyKatException("There is no node called " + value)
    
    @property
    def name(self): return self.__name        

    def __str__(self): return self.name

    def mask(self, n, m, factor):
        id = str(n)+"_"+str(m)
        
        # if the mask is 1 then remove this so it doesn't get 
        # printed as by default the value is 1.0
        if id in self._mask and factor == 1.0:
            del self._mask[id]
                
        self._mask[id] = factor
    
class pd(Detector):

    def __init__(self, name, num_demods, node_name, senstype=None, alternate_beam=False, **kwargs):
        Detector.__init__(self, name, node_name)
        
        self.__num_demods = num_demods
        self.__senstype = senstype
        self.__alternate_beam = alternate_beam
        # create the parameters for all 5 demodulations regardless
        # of how many the user specifies. Later we add properties to
        # those which correspond to the number of demodulations
        
        self.__f1 = Param("f1", self, 0)
        self.__f2 = Param("f2", self, 0)
        self.__f3 = Param("f3", self, 0)
        self.__f4 = Param("f4", self, 0)
        self.__f5 = Param("f5", self, 0)
        
        self.__phi1 = Param("phi1", self, None)
        self.__phi2 = Param("phi2", self, None)
        self.__phi3 = Param("phi3", self, None)
        self.__phi4 = Param("phi4", self, None)
        self.__phi5 = Param("phi5", self, None)
        
        # define new class for assigning new attributes
        cls = type(self)
        self.__class__ = type(cls.__name__, (cls,), {})
    
        self.__set_demod_attrs()
        
    @property
    def senstype(self): return self.__senstype
    @senstype.setter
    def senstype(self,value):
        if value == "": value = None
        
        if value != "S" and value != "N" and value != None: 
            raise pkex.BasePyKatException("Photodiode sensitivity type can either be 'N', 'S' or None.")
            
        self.__senstype = value
        
    @property
    def num_demods(self): return self.__num_demods
    @num_demods.setter
    def num_demods(self, value): 
        if value < 0 or value > 5:
            raise pkex.BasePyKatException("Number of demodulations must be between 0 and 5")
        
        self.__num_demods = value
        self.__set_demod_attrs()
    
    def __get_fphi(self, name):
        return getattr(self, '_'+ self.__class__.__name__ +'__' + name)
    
    def __set_f(self, num, value):
        setattr(self, '_'+ self.__class__.__name__ +'__f' + name, float(value))
    
    def __set_phi(self, num, value):
        if value == None and num != self.num_demods:
            # check if we are setting no phase that this is only on the last
            # demodulation phase.
            raise pkex.BasePyKatException("Only last demodulation phase can be set to None")
        elif isinstance(value, str) and not isinstance(value,float) and value.lower() != "max":
            raise pkex.BasePyKatException("Demodulation phase can only be set to a 'max' or a number (or None if the last demodulation phase)")
            
        setattr(self, '_'+ self.__class__.__name__ +'__phi' + name, value)
        
    def __set_demod_attrs(self):
        """
        For the set number of demodulations the correct number of 
        Parameters are created.
        """
        
        # if there are demodulations present then we want to add
        # the various parameters so they are available for users
        # to play with.
        if self.__num_demods > 0:
            for i in range(1,6):
                name = str(i)
                if i <= self.num_demods:
                    if not hasattr(self, "f"+name):
                        setattr(self.__class__, "f"+name, property(fget=lambda self, i=i: self.__get_fphi('f'+str(i)), fset=lambda self, value, i=i: self.__set_f(str(i), value)))
                    
                    if not hasattr(self, "phi"+name):
                        setattr(self.__class__, "phi"+name, property(fget=lambda self, i=i: self.__get_fphi('phi'+str(i)), fset=lambda self, value, i=i: self.__set_phi(str(i), value)))
                else:
                    if hasattr(self, "f"+name):
                        delattr(self.__class__, "f"+name)
                    if hasattr(self, "phi"+name):
                        delattr(self.__class__, "phi"+name)
        else:
            return
    
    def getFinesseText(self) :
        rtn = []
        
        if self.enabled:
            alt_str = ""
            fphi_str = ""
            
            if self.__alternate_beam:
                alt_str = "*"
                
            for n in range(1, 1+self.num_demods):
                fphi_str += str(self.__getattribute__("f"+str(n)))
                phi_val = self.__getattribute__("phi"+str(n))
                
                if phi_val != None:
                    fphi_str += " " + str(phi_val)
            
            senstype = self.senstype
            
            if senstype == None:
                senstype = ""
                
            rtn.append("pd{0}{1} {2} {3} {4}{5}".format(senstype, self.num_demods, self.name, fphi_str, self.node.name, alt_str))
                
        for p in self._params:
            rtn.extend(p.getFinesseText())
            
        return rtn
            
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

    @property
    def pdtype(self): return self.__pdtype
    @pdtype.setter
    def pdtype(self, value): self.__pdtype = value

    def __init__(self, name, node, senstype="", num_demods=0, demods=[], pdtype=None):        
        Detector.__init__(self, name, node)
        
        if num_demods>2:
            raise NotImplementedError("pd with more than two demodulations not implemented yet")   
            
        self.num_demods = num_demods
        self.senstype = senstype
        self.__pdtype = pdtype


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

            if self.scale != None and self.scale !='':
                rtn.append("scale {1} {0}".format(self.name, self.scale))

            if self.pdtype != None and self.pdtype != '':
                rtn.append("pdtype {0} {1}".format(self.name, self.pdtype))

            if self.noplot:
                rtn.append("noplot {0}".format(self.name))
            
            return rtn
        else:
            return None
    
    def getQGraphicsItem(self):
        if self._svgItem == None:
            self._svgItem = ComponentQGraphicsItem(":/resources/photodiode_red.svg",self,[(-5,11,self.node)])
        
        return self._svgItem    
        
