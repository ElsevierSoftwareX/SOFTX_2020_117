# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 0split()9:09:10 2013

@author: Daniel
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pykat.external.six as six
if six.PY2:
    import exceptions

import abc
from pykat.node_network import *
from pykat.param import Param, AttrParam
from pykat.SIfloat import SIfloat

import pykat.external.six as six
import pykat.exceptions as pkex

import collections
import warnings
import copy

from pykat import USE_GUI, NoGUIException
from pykat.freeze import canFreeze

if USE_GUI:
    import pykat.gui.resources
    from pykat.gui.graphics import *

id_____pykat_class = 0

@canFreeze
class BaseDetector(object) :
    """
    This is a base class for all detectors. Classes Detector1 and Detector2 should be used directly.
    This base class can handled detectors connected to multiple nodes.
    """
    
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
        
    def __init__(self, name, nodes=None, max_nodes=1):

        self._unfreeze()
        
        self.__name = name
        self._svgItem = None
        self._kat = None
        self.noplot = False
        self.enabled = True
        self.tag = None
        self._params = []
        self._mask = {}
        self.__scale = []
        self.__removed = False
        self.noplot = False
        
        self._alternate_beam = []
        self._nodes = []
        self._requested_nodes = []
        
        if nodes != None:
            if isinstance(nodes, (list, tuple)):
                
                if len(nodes) > max_nodes:
                    raise pkex.BasePyKatException("Tried to set too many nodes, %s, maximum number is %i." %(str(nodes),max_nodes))
                    
                for n in nodes:
                    if n[-1]=='*':
                        self._alternate_beam.append(True)
                        n = n[:-1]
                    else:
                        self._alternate_beam.append(False)
                        
                    self._requested_nodes.append(n)
            elif isinstance(nodes, six.string_types):
                # if we don't have a collection
                if nodes[-1]=='*':
                    self._alternate_beam.append(True)
                    nodes = nodes[:-1]
                else:
                    self._alternate_beam.append(False)
                    
                self._requested_nodes.append(nodes)
            else:
                raise pkex.BasePyKatException("Nodes should be a list or tuple of node names or a singular node name as a string.")
    
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
        
        result._freeze()
        return result
                
    def _register_param(self, param):
        self._params.append(param)
        
    def _on_kat_add(self, kat):
        self._kat = kat
        
        for rn in self._requested_nodes:
            if rn != None:
                self._nodes.append(kat.nodes.createNode(rn))
    
    def _on_kat_remove(self):
        self.__removed = True
        
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
        #raise NotImplementedError("This function is not implemented")
        
        if self.noplot:
            return ["noplot " + self.name]
        else:
            return []
        
    def getQGraphicsItem(self):    
        if not USE_GUI:
            raise NoGUIException
            
        return None

    @property
    def removed(self): return self.__removed
    
    @property 
    def scale(self): return self.__scale
    @scale.setter
    def scale(self, value):
        if value not in self.__scale:
            self.__scale = value
    
    @property
    def name(self): return self.__name        

    def __str__(self): return self.name

    def mask(self, n, m, factor):
        _id = str(n)+"_"+str(m)
        
        # if the mask is 1 then remove this so it doesn't get 
        # printed as by default the value is 1.0
        if _id in self._mask and factor == 1.0:
            del self._mask[_id]
                
        self._mask[_id] = factor

    def _set_node(value, index):
        if self._kat is None:
            raise pkex.BasePyKatException("This detector has not been added to a kat object yet")
        else:
            if value[-1] == '*':
                self._alternate_beam[index] = True
                value = value[:-1]
            else:
                self._alternate_beam[index] = False
                
            if value in self._kat.nodes:
                self._nodes[index] = self._kat.nodes[value]
            else:
                raise pkex.BasePyKatException("There is no node called " + value + " in the kat object this detector is attached to.")
    
    def _getScaleCmds(self, rtn):
        if self.scale != None:
            if isinstance(self.scale, six.string_types):
                rtn.append("scale {1} {0}".format(self.name, self.scale))
            elif isinstance(self.scale, (list, tuple)):
                for s in self.scale:
                    rtn.append("scale {1} {0}".format(self.name, s))
            else:
                raise pkex.BasePyKatException("Scale command should either be a list of scales or a single string.")

class Detector0(BaseDetector):
    """
    A detector that attaches to no nodes.
    """
    pass
                        
class Detector1(BaseDetector):
    """
    A detector that attaches to one node.
    """
    @property 
    def node(self): return self._nodes[0]
    @node.setter
    def node(self, value):
        self._set_node(value, 0)      
        
                
class Detector2(BaseDetector):
    """
    A detector that attaches to two node.
    """
    
    @property 
    def node1(self): return self._nodes[0]
    @node1.setter
    def node(self, value):
        self._set_node(value, 0)
        
    @property 
    def node2(self): return self._nodes[1]
    @node2.setter
    def node(self, value):
        self._set_node(value, 1)
                
                
class beam(Detector1):
    
    def __init__(self, name, node_name, frequency=None, alternate_beam=False):
        BaseDetector.__init__(self, name, node_name)
        
        self.alternate_beam = alternate_beam
        self.__f = Param("f", self, frequency)        
    
        self._freeze()
        
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
        
        if len(values) == 3:
            return beam(values[1], node, alternate_beam=alt_beam)
        elif len(values) == 4:
            return beam(values[1], node, alternate_beam=alt_beam, frequency=SIfloat(values[2]))
        else:
            raise pkex.BasePyKatException('Beam detector code "{0}" is not a valid FINESSE command'.format(text))
    
    def getFinesseText(self) :
        rtn = BaseDetector.getFinesseText(self)
        
        if self.alternate_beam:
            alt = '*'
        else:
            alt = ''
        
        if self.f.value is None:
            rtn.append("beam {name} {node}{alt}".format(name=self.name, node=self.node.name, alt=alt))
        else:
            rtn.append("beam {name} {f} {node}{alt}".format(name=self.name, f=str(self.f.value), node=self.node.name, alt=alt))
        
        for p in self._params:
            rtn.extend(p.getFinesseText())
        
        return rtn
        
class cp(Detector0):
    
    def __init__(self, name, cavity, direction, parameter):
        BaseDetector.__init__(self, name, None)
        
        self.cavity = str(cavity)
        self.direction = direction
        self.parameter = parameter
        
        self._freeze()

    @property
    def direction(self): return self.__direction
    @direction.setter
    def direction(self, value):
        if value.lower() not in ["x", "y"]:
            raise pkex.BasePyKatException('Cavity parameter detector direction must be x or y.')
    
        self.__direction = value
        
        
    @property
    def parameter(self): return self.__param
    @parameter.setter
    def parameter(self, value):
        
        params = ["w0","w","zr","z","r","q","finesse","m","stability","loss","length","fsr","fwhm","pole","gouy","fsep","A","B","C","D"]

        if value not in params:
            raise pkex.BasePyKatException('Cavity parameter detector direction must be one of: ' + ", ".join(params))
    
        self.__param = value


    @staticmethod
    def parseFinesseText(text): 
        values = text.split()
        
        if len(values) == 4:
            # For FINESSE version < 2.1
            # Old format the name of the detector is a combination of the arguments
            return cp(values[1] + "_" + values[2] + "_" + values[3], values[1], values[2], values[3])
        elif len(values) == 5:
            return cp(values[1], values[2], values[3], values[4])
        else:
            raise pkex.BasePyKatException('Cavity parameter detector code "{0}" is not a valid FINESSE command'.format(text))
            
    def getFinesseText(self) :
        rtn = BaseDetector.getFinesseText(self)
        
        rtn.append("cp {name} {cavity} {direction} {parameter}".format(name=self.name,
                                                                   cavity=str(self.cavity),
                                                                   direction=self.direction,
                                                                   parameter=self.parameter))
        
        return rtn
        
class xd(Detector0):
    
    def __init__(self, name, component, motion):
        BaseDetector.__init__(self, name, None)
        
        self.component = component
        self.motion = motion
    
        self._freeze()

    @staticmethod
    def parseFinesseText(text): 
        values = text.split()
        
        if len(values) == 4:
            return xd(values[1], values[2], values[3])
        else:
            raise pkex.BasePyKatException('Motion detector code "{0}" is not a valid FINESSE command'.format(text))
            
    def getFinesseText(self) :
        rtn = BaseDetector.getFinesseText(self)
        
        rtn.append("xd {name} {component} {motion}".format(name=self.name,
                                                           component=self.component,
                                                           motion=self.motion))
        
        return rtn
              
        
class ad(Detector1):
    
    def __init__(self, name, frequency, node_name, mode=None, alternate_beam=False):
        BaseDetector.__init__(self, name, node_name)
        self.mode = mode
        self.alternate_beam = alternate_beam
        self.__f = Param("f", self, frequency)
    
        self._freeze()
    
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
        rtn = BaseDetector.getFinesseText(self)
        
        if self.alternate_beam:
            alt = '*'
        else:
            alt = ''
        
        if self.mode is None:
            rtn.append("ad {name} {f} {node}{alt}".format(name=self.name, f=str(self.f.value), node=self.node.name, alt=alt))
        else:
            rtn.append("ad {name} {n} {m} {f} {node}{alt}".format(name=self.name, n=str(self.mode[0]), m=str(self.mode[1]), f=str(self.f.value), node=self.node.name, alt=alt))
            
        for p in self._params:
            rtn.extend(p.getFinesseText())
        
        return rtn

class gouy(Detector1):
    
    def __init__(self, name, direction, spaces):
        BaseDetector.__init__(self, name)
        self.spaces = copy.copy(spaces)
        self.direction = direction
        self.alternate_beam = False
    
        self._freeze()
        
    @property
    def direction(self): return self.__dir
    @direction.setter
    def direction(self, value):
        if value is None or (value != 'x' and value != 'y'):
            raise pkex.BasePyKatException('Direction must be either x or y')
    
        self.__dir = value

    @property
    def spaces(self): return self.__spaces
    @spaces.setter
    def spaces(self, value):

        if value is None or len(value) < 1:
            raise pkex.BasePyKatException('Must be a list of space names')
    
        self.__spaces = value
        
    @staticmethod
    def parseFinesseText(text): 
        values = text.split()

        if len(values) > 3:
            return gouy(str(values[1]), str(values[2]), values[3:])
        else:
            raise pkex.BasePyKatException('Gouy detector code "{0}" is not a valid FINESSE command'.format(text))
            
    def getFinesseText(self) :
        rtn = BaseDetector.getFinesseText(self)

        rtn.append("gouy {name} {dir} {spaces}".format(name=self.name, dir=str(self.direction), spaces=" ".join(self.spaces)))
        
        for p in self._params:
            rtn.extend(p.getFinesseText())
        
        return rtn



class bp(Detector1):
    acceptedParameters = ['w', 'w0', 'z', 'zr', 'g', 'r', 'q', 'Rc']
    
    def __init__(self, name, direction, parameter, node, alternate_beam=False):
        BaseDetector.__init__(self, name, node)
        self.parameter = parameter
        self.direction = direction
        self.alternate_beam = alternate_beam
    
        self._freeze()
        
    @property
    def direction(self): return self.__dir
    @direction.setter
    def direction(self, value):
        if value is None or (value != 'x' and value != 'y'):
            raise pkex.BasePyKatException('Direction must be either x or y')
    
        self.__dir = value

    @property
    def parameter(self): return self.__param
    @parameter.setter
    def parameter(self, value):
        
        if value is None or (value not in self.acceptedParameters) :
            raise pkex.BasePyKatException('Parameter must be one of: %s'%(", ".join(self.acceptedParameters)))
    
        self.__param = value
        
    @staticmethod
    def parseFinesseText(text): 
        values = text.split()
        
        node=values[-1]
        alt_beam = node[-1] == '*'
        
        if len(values) > 3:
            return bp(str(values[1]), str(values[2]), str(values[3]), str(values[4]), alternate_beam=alt_beam)
        else:
            raise pkex.BasePyKatException('Gouy detector code "{0}" is not a valid FINESSE command'.format(text))
            
    def getFinesseText(self) :
        rtn = BaseDetector.getFinesseText(self)

        if self.alternate_beam:
            alt = "*"
        else:
            alt = ""
            
        rtn.append("bp {name} {dir} {param} {node}{alt}".format(name=self.name, dir=str(self.direction), param=self.parameter, node=self.node.name, alt=alt))
        
        for p in self._params:
            rtn.extend(p.getFinesseText())
        
        return rtn

id___2 = 0

class pd(Detector1):

    def __new__(cls, *args, **kwargs):
        # This creates an instance specific class for the component
        # this enables us to add properties to instances rather than
        # all classes
        global id___2
        id___2 += 1
        cnew_name = str("%s.%s_%i" % (cls.__module__, cls.__name__, id___2))
    
        cnew = type(cnew_name, (cls,), {})
    
        return object.__new__(cnew)

    def __deepcopy__(self, memo):
        """
        When deep copying a kat object we need to take into account
        the instance specific properties.
        """
        
        # Here we create a copy of this object based of the base class
        # of this one, otherwise we're making a copy of a copy of a copy...
        result = self.__class__.__new__(self.__class__.__base__)
        result.__dict__ = copy.deepcopy(self.__dict__, memo)
        
        result.__set_demod_attrs()
        
        return result

    def __init__(self, name=None, num_demods=1, node_name=None, senstype=None, alternate_beam=False, pdtype=None, **kwargs):
        BaseDetector.__init__(self, name, node_name)
        
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
        
        self.__phase1 = Param("phase1", self, None)
        self.__phase2 = Param("phase2", self, None)
        self.__phase3 = Param("phase3", self, None)
        self.__phase4 = Param("phase4", self, None)
        self.__phase5 = Param("phase5", self, None)
        
        fs = [self.__f1, self.__f2, self.__f3, self.__f4, self.__f5]
        ps = [self.__phase1, self.__phase2, self.__phase3, self.__phase4, self.__phase5]

        for i in range(num_demods):
            f = 'f{0}'.format(i+1)
            #print("i {0} fs {1} f {2} keys {3}".format(i,len(fs),f, kwargs.keys()))
            
            if f in kwargs:
                fs[i].value = kwargs[f]
            else:
                raise pkex.BasePyKatException("Missing demodulation frequency {0} (f{0})".format(i+1))    
        
            p = 'phase{0}'.format(i+1)
            
            if p in kwargs:
                if kwargs[p] is None and i<num_demods-1:
                    raise pkex.BasePyKatException("Missing demodulation phase {0} (phase{0})".format(i+1))
                    
                ps[i].value = kwargs[p]
            elif i<num_demods-1:
                raise pkex.BasePyKatException("Missing demodulation phase {0} (phase{0})".format(i+1))
   
        self.__set_demod_attrs()
    
        self._freeze()
                
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
        value = SIfloat(value)
        
        p = getattr(self, '_pd__f' + num)
        p.value = value
    
    def __set_phi(self, num, value):
        value = SIfloat(value)
        
        if value is None and num != self.num_demods:
            # check if we are setting no phase that this is only on the last
            # demodulation phase.
            raise pkex.BasePyKatException("Only last demodulation phase can be set to None")
        elif isinstance(value, six.string_types) and not isinstance(value,float) and value.lower() != "max":
            raise pkex.BasePyKatException("Demodulation phase can only be set to a 'max' or a number (or None if the last demodulation phase)")
          
        p = getattr(self, '_pd__phase' + num)  
        p.value = value
        
    def __set_demod_attrs(self):
        """
        For the set number of demodulations the correct number of 
        Parameters are created.
        """
        self._unfreeze()
        
        # if there are demodulations present then we want to add
        # the various parameters so they are available for users
        # to play with.
        if self.__num_demods > 0:
            for i in range(1,6):
                name = str(i)
                
                if i <= self.num_demods:
                    if not hasattr(self, "f"+name):
                        fget = lambda self, i=i:        self.__get_fphi('f'+str(i))
                        fset = lambda self, value, i=i: self.__set_f(str(i), value)
                        
                        setattr(self.__class__, "f"+name, property(fget=fget, fset=fset))
                    
                    if not hasattr(self, "phase"+name):
                        setattr(self.__class__, "phase"+name, property(fget=lambda self, i=i: self.__get_fphi('phase'+str(i)), fset=lambda self, value, i=i: self.__set_phi(str(i), value)))
                else:
                    if hasattr(self, "f"+name):
                        delattr(self.__class__, "f"+name)
                        
                    if hasattr(self, "phase"+name):
                        delattr(self.__class__, "phase"+name)
        
        self._freeze()
    
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
            dict['f{0}'.format(i+1)] = SIfloat(f[i])
        for i in range(len(phs)):
            if phs[i] == "max":
                dict['phase{0}'.format(i+1)] = "max"
            else:
                dict['phase{0}'.format(i+1)] = SIfloat(phs[i])
            
        node = values[-1]
        alt_beam = node[-1] == '*'
        
        if alt_beam:
            node = node[0:-1]
        
        return pd(values[1], demods, node, senstype=senstype, alternate_beam=alt_beam, **dict)

        
    def getFinesseText(self) :
        rtn = BaseDetector.getFinesseText(self)
        
        if self.enabled:
            alt_str = ""
            fphi_str = ""
            
            if self.alternate_beam:
                alt_str = "*"
                
            for n in range(1, 1+self.num_demods):
                _f = self.__getattribute__("f"+str(n))
                
                if _f == "$fs":
                    fphi_str += " $fs"
                else:
                    fphi_str += " {0:.16g}".format(float(_f))
                    
                phi_val = self.__getattribute__("phase"+str(n))
                
                if phi_val != None:
                    if type(phi_val) == float:
                        fphi_str += " {0:.16g}".format(float(phi_val))
                    else:
                        fphi_str += " {0}".format(phi_val)
            
            senstype = self.senstype
            
            if senstype is None:
                senstype = ""
            
            rtn.append("pd{0}{1} {2}{3} {4}{5}".format(senstype, self.num_demods, self.name, fphi_str, self.node.name, alt_str))

            self._getScaleCmds(rtn)
            
            if self.pdtype != None:
                rtn.append("pdtype {0} {1}".format(self.name, self.pdtype))
                
            for p in self._params:
                rtn.extend(p.getFinesseText())
            
        return rtn
  
class qnoised(pd):
    
    def __init__(self, name, num_demods, node_name, alternate_beam=False, pdtype=None, **kwargs):
        super(qnoised, self).__init__(name, num_demods, node_name, alternate_beam=alternate_beam, pdtype=pdtype, **kwargs)
    
        self._unfreeze()
        
        self.__homangle = AttrParam("homangle", self, None)
        
        self._freeze()
    
    @property
    def homangle(self): return self.__homangle
    @homangle.setter
    def homangle(self, value): self.__homangle.value = value
    
    @pd.pdtype.setter
    def pdtype(self, value):
        raise pkex.BasePyKatException("Setting pdtype is not possible with qnoised detectors")
    
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
            dict['phase{0}'.format(i+1)] = phs[i]
            
        node = values[-1]
        alt_beam = node[-1] == '*'
        
        if alt_beam:
            node = node[0:-1]
        
        if values[0].endswith('S'):
            sens='S'
        elif values[0].endswith('N'):
            sens='N'
        else:
            sens=None
            
        return qnoised(values[1], demods, node, senstype=sens, alternate_beam=alt_beam, **dict)
    
    def getFinesseText(self) :
        rtn = BaseDetector.getFinesseText(self)
        
        if self.enabled:
            alt_str = ""
            fphi_str = ""
            
            if self.alternate_beam:
                alt_str = "*"
            
            for n in range(1, 1+self.num_demods):
                _f = self.__getattribute__("f"+str(n))
                
                if _f == "$fs":
                    fphi_str += " $fs"
                else:
                    fphi_str += " {0:.16g}".format(float(_f))
                    
                phi_val = self.__getattribute__("phase"+str(n))
                
                if phi_val != None:
                    if type(phi_val) == float:
                        fphi_str += " {0:.16g}".format(float(phi_val))
                    else:
                        fphi_str += " " + str(phi_val)
            
            senstype = self.senstype
            
            if senstype is None:
                senstype = ""
                
            rtn.append("qnoised{5} {0} {1} {2} {3}{4}".format(self.name, self.num_demods, fphi_str, self.node.name, alt_str, senstype))

            self._getScaleCmds(rtn)
                
            for p in self._params:
                rtn.extend(p.getFinesseText())
            
        return rtn

class qshot(pd):
    
    def __init__(self, name, num_demods, node_name, alternate_beam=False, **kwargs):
        super(qshot, self).__init__(name, num_demods, node_name, alternate_beam=alternate_beam, pdtype=None, **kwargs)     

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
            dict['phase{0}'.format(i+1)] = phs[i]
            
        node = values[-1]
        alt_beam = node[-1] == '*'
        
        if alt_beam:
            node = node[0:-1]
            
        if values[0].endswith('S'):
            sens='S'
        elif values[0].endswith('N'):
            sens='N'
        else:
            sens=None
        
        return qshot(values[1], demods, node, senstype=sens, alternate_beam=alt_beam, **dict)
    
    def getFinesseText(self) :
        rtn = BaseDetector.getFinesseText(self)
        
        if self.enabled:
            alt_str = ""
            fphi_str = ""
            
            if self.alternate_beam:
                alt_str = "*"
                
            for n in range(1, 1+self.num_demods):
                _f = self.__getattribute__("f"+str(n))
                
                if _f == "$fs":
                    fphi_str += " $fs"
                else:
                    fphi_str += " {0:.16g}".format(float(_f))
                    
                phi_val = self.__getattribute__("phase"+str(n))
                
                if phi_val != None:
                    if type(phi_val) == float:
                        fphi_str += " {0:.16g}".format(float(phi_val))
                    else:
                        fphi_str += " " + str(phi_val)
            
            senstype = self.senstype
            
            if senstype is None:
                senstype = ""
                
            rtn.append("qshot{5} {0} {1} {2} {3}{4}".format(self.name, self.num_demods, fphi_str, self.node.name, alt_str,senstype))

            self._getScaleCmds(rtn)
                
            for p in self._params:
                rtn.extend(p.getFinesseText())
            
        return rtn
    
    
class hd(Detector2):
    
    def __init__(self, name, phase, node1_name, node2_name):
        BaseDetector.__init__(self, name, (node1_name, node2_name), max_nodes=2)
    
        self.__phase = Param("phase", self, phase)

        self._freeze()
        
    @property
    def phase(self): return self.__phase
    @phase.setter
    def phase(self, value): self.__phase.value = value
    
    def parseAttributes(self, values):
        raise pkex.BasePyKatException("hd detector %s has no attributes to set" % self.name)
    
    @staticmethod
    def parseFinesseText(text): 
        values = text.split()
        
        return hd(values[1], float(values[2]), str(values[3]), str(values[4]))
    
    def getFinesseText(self):
        rtn = BaseDetector.getFinesseText(self)
        
        if self.enabled:   
            n1 = self.node1.name
            n2 = self.node2.name
            
            if self._alternate_beam[0]: n1 += "*"
            if self._alternate_beam[1]: n2 += "*"
            
            rtn.append("hd {0} {1} {2} {3}".format(self.name, self.phase, n1, n2))

            self._getScaleCmds(rtn)
                
            for p in self._params:
                rtn.extend(p.getFinesseText())
            
        return rtn
        
class qhd(Detector2):
    
    def __init__(self, name, phase, node1_name, node2_name, sensitivity=""):
        BaseDetector.__init__(self, name, (node1_name, node2_name), max_nodes=2)
    
        self.__phase = Param("phase", self, phase)
        self.sensitivity = sensitivity
    
        self._freeze()
        
    @property
    def phase(self): return self.__phase
    @phase.setter
    def phase(self, value): self.__phase.value = value
    
    @property
    def sensitivity(self): 
        return self.__sensitivity
    @sensitivity.setter
    def sensitivity(self, value):
        if value == 'S' or value == 'N':
            self.__sensitivity = value
        elif value is None or value == '':
            self.__sensitivity = ""
        else:
            raise pkex.BasePyKatException("qhd (%s) sensitivity option '%s' is not available, use either 'S' or 'N'." % (self.name, value))
        
    
    def parseAttributes(self, values):
        raise pkex.BasePyKatException("hd detector %s has no attributes to set" % self.name)
    
    @staticmethod
    def parseFinesseText(text): 
        values = text.split()
        
        sens = values[0][-1]
        
        if sens != 'S' and sens != 'N':
            sens = None
            
        return qhd(values[1], float(values[2]), str(values[3]), str(values[4]), sensitivity = sens)
    
    def getFinesseText(self):
        rtn = BaseDetector.getFinesseText(self)
        
        if self.enabled:   
            n1 = self.node1.name
            n2 = self.node2.name
            
            if self._alternate_beam[0]: n1 += "*"
            if self._alternate_beam[1]: n2 += "*"
            
            rtn.append("qhd{4} {0} {1} {2} {3}".format(self.name, self.phase, n1, n2, self.sensitivity))

            self._getScaleCmds(rtn)
                
            for p in self._params:
                rtn.extend(p.getFinesseText())
            
        return rtn
