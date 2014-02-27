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
from pykat.param import Param, putter
import pykat.exceptions as pkex
from collections import namedtuple

class Command(object):
    def __init__(self):
        self.tag = None

    def getFinesseText(self):
        """ Base class for individual finesse optical components """
        raise NotImplementedError("This function is not implemented")

    @staticmethod
    def parseFinesseText(text):
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

class gauss(object):
    @staticmethod
    def parseFinesseText(text, kat):
        values = text.split()

        if not values[0].startswith("gauss") or (len(values) != 6 and len(values != 8)):
            raise exceptions.RuntimeError("'{0}' not a valid Finesse gauss command".format(text))        
        
        name = values[1]
        component = values[2]
        node = values[3]
        
        
        
        if values[0].endswith("**"):
            if len(values) == 6:
                print ""
            elif len(values) == 8:
                print ""
        elif values[0].endswith("*"):
            if len(values) == 6:
                print ""
            elif len(values) == 8:
                print ""
        else:
            if len(values) == 6:
                print ""
            elif len(values) == 8:
                print ""
            
class tf(Command):
    fQ = namedtuple('fQ', 'f', 'Q')
    
    def __init__(self, name, poles, zeros):
        pass
      
class xaxis(Command):
    """
    The xaxis object is a unique object to each pykat.finesse.kat instance. It provides
    and interface to the xaxis command in FINESSE.
    """
    
    def __init__(self, scale, limits, param, steps, comp=None, axis_type="xaxis"):
        """
        Typical usage:
            xaxis(["lin" or "log"], [upper, lower], param, steps)
            
        param must be an object of the type pykat.param.Param whose
        isPutable() member returns true.
        
        steps is the number of points to compute between upper and lower limits.
        """
        self._axis_type = axis_type

        self.x = putter("x1")
        self.mx = putter("mx1")

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

        self.limits = numpy.array(SIfloat(limits)).astype(float)

        if steps <= 0 :
            raise exceptions.ValueError("steps value should be > 0")

        self.steps = int(steps)

        if isinstance(param, str):
            self.__param = param
            if comp == None:
                raise pkex.BasePyKatException("If parameter is set with a string, the comp argument must set the component name")
                
            self.__comp = comp
        elif not isinstance(param, Param) :
            raise pkex.BasePyKatException("param argument is not of type Param")
        else:
            self.__param = param
            self.__comp = param._owner.name

    @property
    def param(self): return self.__param
    @param.setter
    def param(self, value):
	if not isinstance(value, Param):
		raise pkex.BasePyKatException("param argument is not of type Param")
	else:
		self.__param = value
		self.__comp = value._owner.name

    @staticmethod
    def parseFinesseText(text):
        values = text.split()

        if values[0] != "xaxis" and values[0] != "xaxis*":
            raise exceptions.RuntimeError("'{0}' not a valid Finesse xaxis command".format(text))

        axis_type = values[0]

        values.pop(0) # remove initial value

        if len(values) != 6:
            raise exceptions.RuntimeError("xaxis Finesse code format incorrect '{0}'".format(text))

        return xaxis(values[2], [values[3], values[4]], values[1], values[5], comp=values[0], axis_type=axis_type)

    def getFinesseText(self):
        # store either the component name of the string provided
        comp_name = self.__comp.name if isinstance(self.__comp, Component) else self.__comp
        param_name = self.__param.name if isinstance(self.__param, Param) else self.__param

        return '{axis_type} {0} {1} {2} {3} {4} {5}'.format(
                comp_name, param_name, self.scale,
                min(self.limits), max(self.limits), self.steps, axis_type=self._axis_type);

class x2axis(xaxis):
    def __init__(self, scale, limits, param, steps, comp=None):
        xaxis.__init__(self, scale, limits, param, steps, comp=comp, axis_type="x2axis")
        self.x = putter("x2")
        self.mx = putter("mx2")

    @staticmethod
    def parseFinesseText(text):
        values = text.split()

        if values[0] != "x2axis" and values[0] != "x2axis*":
            raise exceptions.RuntimeError("'{0}' not a valid Finesse xaxis command".format(text))

        axis_type = values[0]

        values.pop(0) # remove initial value

        if len(values) != 6:
            raise exceptions.RuntimeError("xaxis Finesse code format incorrect '{0}'".format(text))

        return x2axis(values[2], [values[3], values[4]], values[1], values[5], comp=values[0])
