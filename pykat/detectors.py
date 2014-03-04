# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 0split()9:09:10 2013

@author: Daniel
"""
import exceptions
import pykat.gui.resources

from pykat.utils import *
from pykat.gui.graphics import *
from pykat.node_network import *
from pykat.param import Param

import pykat.exceptions as pkex
import warnings

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
        self.__scale = None
        self.__removed = False
        
        if node != None:
            if node[-1]=='*':
                self._alternate_beam = True
                node=node[:-1]
            
            self.__requested_node = node
    
    def _register_param(self, param):
        self._params.append(param)
        
    def _on_kat_add(self, kat):
        self._kat = kat
        
        if self.__requested_node != None:
            self.__node = kat.nodes.createNode(self.__requested_node)
    
    def remove(self):
        if self.__removed:
            raise pkex.BasePyKatException("{0} has already been marked as removed".format(self.name))
        else:
            self._kat.remove(self)
    
        self.__removed = True
        
    @staticmethod
    def parseFinesseText(text):    
        raise NotImplementedError("This function is not implemented")
        
    def getFinesseText(self):
        """ Base class for individual finesse optical components """    
        raise NotImplementedError("This function is not implemented")
        
    def getQGraphicsItem(self):    
        return None

        
    @property
    def removed(self): return self.__removed
    
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

class ad(Detector):
    
    def __init__(self, name, frequency, node_name, mode=None, alternate_beam=False):
        Detector.__init__(self, name, node_name)
        self.mode = mode
        self.alternate_beam = alternate_beam
        self.__f = Param("f", self, frequency)
    
    @property
    def mode(self): return self.__mode
    @mode.setter
    def mode(self, value):
        if value != None and len(value) != 2:
            raise pkex.BasePyKatException('Mode must be a container of length 2, first element the x mode and second the y mode')
    
        self.__mode = value
        
    @property
    def f(self): return self.__f
    @f.setter
    def f(self, value): 
        self.__f.value = value
        
    @staticmethod
    def parseFinesseText(text): 
        values = text.split()

        node=values[-1]
        alt_beam = node[-1] == '*'
        if len(values) == 6:
            return ad(values[1], values[4], values[5], mode = [int(values[2]), int(values[3])], alternate_beam=alt_beam)
        elif len(values) == 4:
            return ad(values[1], values[2], values[3], alternate_beam=alt_beam)
        else:
            raise pkex.BasePyKatException('Amplitude detector code "{0}" is not a valid FINESSE command'.format(text))
            
    def getFinesseText(self) :
        rtn = []
        
        if self.alternate_beam:
            alt = '*'
        else:
            alt = ''
        
        if self.mode == None:
            rtn.append("ad {name} {f} {node}{alt}".format(name=self.name, f=str(self.f.value), node=self.node.name, alt=alt))
        else:
            rtn.append("ad {name} {n} {m} {f} {node}{alt}".fomat(name=self.name, n=str(self.mode[0]), m=str(self.mode[1]), f=str(self.f.value), node=self.node.name, alt=alt))
            
        for p in self._params:
            rtn.extend(p.getFinesseText())
        
        return rtn
        
class pd(Detector):

    def __init__(self, name, num_demods, node_name, senstype=None, alternate_beam=False, pdtype=None, **kwargs):
        Detector.__init__(self, name, node_name)
        
        self.__num_demods = num_demods
        self.__senstype = senstype
        self.alternate_beam = alternate_beam
        self.__pdtype = pdtype

        # create the parameters for all 5 demodulations regardless
        # of how many the user specifies. Later we add properties to
        # those which correspond to the number of demodulations
        
        self.__f1 = Param("f1", self, None)
        self.__f2 = Param("f2", self, None)
        self.__f3 = Param("f3", self, None)
        self.__f4 = Param("f4", self, None)
        self.__f5 = Param("f5", self, None)
        
        self.__phi1 = Param("phi1", self, None)
        self.__phi2 = Param("phi2", self, None)
        self.__phi3 = Param("phi3", self, None)
        self.__phi4 = Param("phi4", self, None)
        self.__phi5 = Param("phi5", self, None)
        
        fs = [self.__f1, self.__f2, self.__f3, self.__f4, self.__f5]
        ps = [self.__phi1, self.__phi2, self.__phi3, self.__phi4, self.__phi5]
        
        for i in range(num_demods):
            f = 'f{0}'.format(i+1)
            
            if f in kwargs:
                fs[i].value = kwargs[f]
            else:
                raise pkex.BasePyKatException("Missing demodulation frequency {0} (f{0})".format(i+1))    
        
            p = 'phi{0}'.format(i+1)
            
            if p in kwargs:
                if kwargs[p] == None and i<num_demods-1:
                    raise pkex.BasePyKatException("Missing demodulation phase {0} (phi{0})".format(i+1))
                    
                ps[i].value = kwargs[p]
            elif i<num_demods-1:
                raise pkex.BasePyKatException("Missing demodulation phase {0} (phi{0})".format(i+1))
        
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

    @property
    def pdtype(self): return self.__pdtype
    @pdtype.setter
    def pdtype(self, value): self.__pdtype = value
    
    def __get_fphi(self, name):
        return getattr(self, '_pd__' + name)
    
    def __set_f(self, num, value):
        setattr(self, '_pd__f' + num, float(value))
    
    def __set_phi(self, num, value):
        if value == None and num != self.num_demods:
            # check if we are setting no phase that this is only on the last
            # demodulation phase.
            raise pkex.BasePyKatException("Only last demodulation phase can be set to None")
        elif isinstance(value, str) and not isinstance(value,float) and value.lower() != "max":
            raise pkex.BasePyKatException("Demodulation phase can only be set to a 'max' or a number (or None if the last demodulation phase)")
            
        setattr(self, '_'+ self.__class__.__name__ +'__phi' + num, value)
        
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
    
    @staticmethod
    def parseFinesseText(text): 
        values = text.split()
        demods = 0
        senstype = None

        if len(values[0]) == 4:
            senstype = values[0][2]
            demods = int(values[0][3])
        elif len(values[0]) == 3:
            demods = int(values[0][2])
        elif len(values[0]) != 2:
            raise pkex.BasePyKatException("Photodiode code format incorrect '{0}' (1)".format(text))
        
        if len(values) <= 3 and demods > 0:
            raise pkex.BasePyKatException("Photodiode code format incorrect '{0}' (2)".format(text))
        elif len(values) > 3 and demods == 0:
            raise pkex.BasePyKatException("Photodiode code format incorrect '{0}' (3)".format(text))
            
        num_f_phs = len(values) - 3
        expected_f_phs = demods * 2
        
        if not (num_f_phs == expected_f_phs or num_f_phs == expected_f_phs-1):
            raise pkex.BasePyKatException("Photodiode code format incorrect '{0}' (4)".format(text))
        
        f = values[2:len(values)-1:2]    
        phs = values[3:len(values)-1:2]
        
        dict = {}
        
        for i in range(len(f)):
            dict['f{0}'.format(i+1)] = f[i]
        for i in range(len(phs)):
            dict['phi{0}'.format(i+1)] = phs[i]
            
        node = values[-1]
        alt_beam = node[-1] == '*'
        
        if alt_beam:
            node = node[0:-1]
        
        return pd(values[1], demods, node, senstype=senstype, alternate_beam=alt_beam, **dict)

        
    def getFinesseText(self) :
        rtn = []
        
        if self.enabled:
            alt_str = ""
            fphi_str = ""
            
            if self.alternate_beam:
                alt_str = "*"
                
            for n in range(1, 1+self.num_demods):
                fphi_str += " " + str(self.__getattribute__("f"+str(n)))
                phi_val = self.__getattribute__("phi"+str(n))
                
                if phi_val != None:
                    fphi_str += " " + str(phi_val)
            
            senstype = self.senstype
            
            if senstype == None:
                senstype = ""
            
            rtn.append("pd{0}{1} {2}{3} {4}{5}".format(senstype, self.num_demods, self.name, fphi_str, self.node.name, alt_str))

            if self.scale != None:
                rtn.append("scale {1} {0}".format(self.name, self.scale))

            if self.pdtype != None:
                rtn.append("pdtype {0} {1}".format(self.name, self.pdtype))
                
            for p in self._params:
                rtn.extend(p.getFinesseText())
            
        return rtn
  
class qnoised(pd):
    
    def __init__(self, name, num_demods, node_name, alternate_beam=False, pdtype=None, **kwargs):
        super(qnoised, self).__init__(name, num_demods, node_name, alternate_beam=alternate_beam, pdtype=pdtype, **kwargs)
    
        self.__homangle = AttrParam("homangle", self, None)
    
    @property
    def homangle(self): return self.__homangle
    @homangle.setter
    def homangle(self, value): self.__homangle.value = value
    
    @pd.pdtype.setter
    def pdtype(self, value):
        raise pkex.BasePyKatException("Setting pdtype is not possible with qnoised detectors")
    
    @pd.senstype.setter
    def senstype(self,value):
        raise pkex.BasePyKatException("qnoised detector has no sensitvity type")
    
    def parseAttributes(self, values):
        
        for key in values.keys():
            if key in ["homangle"]:
                self.__homangle.value = values[key]
            else:
                raise pkex.BasePyKatException("No attribute {0} for qnoised".format(key))
    
    @staticmethod
    def parseFinesseText(text): 
        values = text.split()

        if len(values) <= 3:
            raise pkex.BasePyKatException("qnoised code format incorrect '{0}' (2)".format(text))
            
        demods = int(values[2])
        
        if len(values) <= 4 and demods > 0:
            raise pkex.BasePyKatException("qnoised code format incorrect '{0}' (2)".format(text))
        elif len(values) > 4 and demods == 0:
            raise pkex.BasePyKatException("qnoised code format incorrect '{0}' (3)".format(text))
            
        num_f_phs = len(values) - 4
        expected_f_phs = demods * 2
        
        if not (num_f_phs == expected_f_phs or num_f_phs == (expected_f_phs-1)):
            raise pkex.BasePyKatException("qnoised code format incorrect '{0}' (4)".format(text))
        
        f = values[3:len(values)-1:2]    
        phs = values[4:len(values)-1:2]
        
        dict = {}
        
        for i in range(len(f)):
            dict['f{0}'.format(i+1)] = f[i]
        for i in range(len(phs)):
            dict['phi{0}'.format(i+1)] = phs[i]
            
        node = values[-1]
        alt_beam = node[-1] == '*'
        
        if alt_beam:
            node = node[0:-1]
        
        return qnoised(values[1], demods, node, alternate_beam=alt_beam, **dict)
    
    def getFinesseText(self) :
        rtn = []
        
        if self.enabled:
            alt_str = ""
            fphi_str = ""
            
            if self.alternate_beam:
                alt_str = "*"
                
            for n in range(1, 1+self.num_demods):
                fphi_str += " " + str(self.__getattribute__("f"+str(n)))
                phi_val = self.__getattribute__("phi"+str(n))
                
                if phi_val != None:
                    fphi_str += " " + str(phi_val)
            
            senstype = self.senstype
            
            if senstype == None:
                senstype = ""
                
            rtn.append("qnoised {0} {1} {2} {3}{4}".format(self.name, self.num_demods, fphi_str, self.node.name, alt_str))

            if self.scale != None:
                rtn.append("scale {1} {0}".format(self.name, self.scale))
                
            for p in self._params:
                rtn.extend(p.getFinesseText())
            
        return rtn

class qnoised(pd):
    
    def __init__(self, name, num_demods, node_name, alternate_beam=False, **kwargs):
        super(qnoised, self).__init__(name, num_demods, node_name, alternate_beam=alternate_beam, pdtype=None, senstype=None, **kwargs)
    
    @pd.pdtype.setter
    def pdtype(self, value):
        raise pkex.BasePyKatException("Setting pdtype is not possible with qshot detectors")
    
    @pd.senstype.setter
    def senstype(self,value):
        raise pkex.BasePyKatException("qshot detector has no sensitvity type")
    
    @staticmethod
    def parseFinesseText(text): 
        values = text.split()

        if len(values) <= 3:
            raise pkex.BasePyKatException("qshot code format incorrect '{0}' (2)".format(text))
            
        demods = int(values[2])
        
        if len(values) <= 4 and demods > 0:
            raise pkex.BasePyKatException("qshot code format incorrect '{0}' (2)".format(text))
        elif len(values) > 4 and demods == 0:
            raise pkex.BasePyKatException("qshot code format incorrect '{0}' (3)".format(text))
            
        num_f_phs = len(values) - 4
        expected_f_phs = demods * 2
        
        if not (num_f_phs == expected_f_phs or num_f_phs == (expected_f_phs-1)):
            raise pkex.BasePyKatException("qshot code format incorrect '{0}' (4)".format(text))
        
        f = values[3:len(values)-1:2]    
        phs = values[4:len(values)-1:2]
        
        dict = {}
        
        for i in range(len(f)):
            dict['f{0}'.format(i+1)] = f[i]
        for i in range(len(phs)):
            dict['phi{0}'.format(i+1)] = phs[i]
            
        node = values[-1]
        alt_beam = node[-1] == '*'
        
        if alt_beam:
            node = node[0:-1]
        
        return qnoised(values[1], demods, node, alternate_beam=alt_beam, **dict)
    
    def getFinesseText(self) :
        rtn = []
        
        if self.enabled:
            alt_str = ""
            fphi_str = ""
            
            if self.alternate_beam:
                alt_str = "*"
                
            for n in range(1, 1+self.num_demods):
                fphi_str += " " + str(self.__getattribute__("f"+str(n)))
                phi_val = self.__getattribute__("phi"+str(n))
                
                if phi_val != None:
                    fphi_str += " " + str(phi_val)
            
            senstype = self.senstype
            
            if senstype == None:
                senstype = ""
                
            rtn.append("qshot {0} {1} {2} {3}{4}".format(self.name, self.num_demods, fphi_str, self.node.name, alt_str))

            if self.scale != None:
                rtn.append("scale {1} {0}".format(self.name, self.scale))
                
            for p in self._params:
                rtn.extend(p.getFinesseText())
            
        return rtn
        
def xd(Detector):

    def __init__(self, name, node_name, component, motion):
        Detector.__init__(name, None)
        
        self.__motion = motion
        self.__component = component
        
    @property
    def motion(self): return self.__motion
    
    @property
    def component(self): return self.__component
    
    @staticmethod
    def parseFinesseText(text): 
        values = text.split()

        if len(values) != 4:
            raise pkex.BasePyKatException("Motion detector command format incorrect '{0}' (2)".format(text))
            
        return xd(values[1], values[2], values[3])
    
    def getFinesseText(self) :
        rtn = []
        
        if self.enabled:
            rtn.append("xd {0} {1} {2}".format(self.name, self.component, self.motion))

            if self.scale != None:
                rtn.append("scale {1} {0}".format(self.name, self.scale))
                
            for p in self._params:
                rtn.extend(p.getFinesseText())
            
        return rtn
    
    