# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 11:58:09 2013

@author: Daniel
"""
import numpy
from numpy import min,max
import exceptions
from components import *
from structs import *

class Command(object):
    def getFinesseText(self):
        """ Base class for individual finesse optical components """    
        raise NotImplementedError("This function is not implemented")
    
    def _on_kat_add(self, kat):
        """
        Called when this component has been added to a kat object
        """
        self._kat = kat
        
class cavity(Command):
    def __init__(self, name, c1, n1, c2, n2):
        self.__name = name
        self.__c1 = c1
        self.__c2 = c2
        self.__n1 = n1
        self.__n2 = n2
        
    def getFinesseText(self):
        return 'cav {0} {1} {2} {3} {4}'.format(self.__name, self.__c1, self.__n1, self.__c2, self.__n2);
    
    
class xaxis(Command):
    
    def __init__(self, scale, limits, comp, param, steps):
        
        if scale == "lin":
            scale = Scale.linear
        elif scale == "log":
            scale = Scale.logarithmic
        elif isinstance(scale, str): 
            # else we have a string but not a recognisable one
            raise exceptions.ValueError("scale argument '{0}' is not valid, must be 'lin' or 'log'".format(scale))
         
        if scale != Scale.linear and scale != Scale.logarithmic:
            raise exceptions.ValueError("scale is not Scale.linear or Scale.logarithmic")            
            
        self.scale = scale
        
        if numpy.size(limits) != 2 :
            raise exceptions.ValueError("limits input should be a 2x1 vector of limits for the xaxis")
            
        self.limits = limits
        
        if steps <= 0 :
            raise exceptions.ValueError("steps value should be > 0")            
            
        self.steps = steps
        
        # if entered component is a string then just store and use that
        if isinstance(comp, str):
            self.__comp = comp
        elif not isinstance(comp, Component):
            raise exceptions.ValueError("{0} has not been added".format(comp))
        else:
            self.__comp = comp
        
        if isinstance(param, str):
            self.__param = param
        elif not isinstance(param, Param) :
            raise exceptions.ValueError("param argument is not of type Param")
        else:
            self.__param = param
        
    def getFinesseText(self):
        # store either the component name of the string provided
        comp_name = self.__comp.name if isinstance(self.__comp, Component) else self.__comp
        param_name = self.__param.name if isinstance(self.__param, Param) else self.__param
        
        return 'xaxis {0} {1} {2} {3} {4} {5}'.format(
                comp_name, param_name, self.scale,
                min(self.limits), max(self.limits), self.steps);
                
        