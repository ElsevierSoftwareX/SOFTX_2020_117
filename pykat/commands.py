# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 11:58:09 2013

@author: Daniel
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy
from numpy import min,max
import pykat.external.six as six
if six.PY2:
	import exceptions
from pykat.components import *
from pykat.structs import *

from pykat.param import Param, putter
import pykat.exceptions as pkex
from collections import namedtuple
from pykat.optics.gaussian_beams import beam_param

class Command(object):
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, name, unique):
        self.__unique = unique
        self.tag = None
        self.__removed = False
        self.__name = name.strip("*")
        
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

    def remove(self):
        self._kat.remove(self)
        self.__removed = True
    
    @property
    def name(self): return self.__name
    
    @property
    def removed(self): return self.__removed
    
class cavity(Command):
    def __init__(self, name, c1, n1, c2, n2):
        Command.__init__(self, name, False)
        
        self.__c1 = c1
        self.__c2 = c2
        self.__n1 = n1
        self.__n2 = n2
        
        self.enabled = True

    def getFinesseText(self):
        if self.enabled:
            return 'cav {0} {1} {2} {3} {4}'.format(self.name, self.__c1.name, self.__n1.name, self.__c2.name, self.__n2.name);
        else:
            return None

    @staticmethod
    def parseFinesseText(line, kat):
        v = line.split()
        
        if len(v) != 6:
            raise pkex.BasePyKatException("cav command format `{0}` is incorrect".format(line))
        
        if v[2] not in kat.components:
            raise pkex.BasePyKatException("cav command `{0}` refers to component `{1}` which does not exist".format(line, v[2]))
        
        if v[4] not in kat.components:
            raise pkex.BasePyKatException("cav command `{0}` refers to component `{1}` which does not exist".format(line, v[4]))
        
        if v[3] not in kat.nodes.getNodes():
            raise pkex.BasePyKatException("cav command `{0}` refers to node `{1}` which does not exist".format(line, v[3]))
        
        if v[5] not in kat.nodes.getNodes():
            raise pkex.BasePyKatException("cav command `{0}` refers to node `{1}` which does not exist".format(line, v[5]))
        
        c1 = getattr(kat, v[2])
        c2 = getattr(kat, v[4])
        
        n1 = getattr(kat.nodes, v[3])
        n2 = getattr(kat.nodes, v[5])
        
        if not hasattr(c1, n1.name):
            raise pkex.BasePyKatException("cav command `{0}`: node `{1}` is not attached to `{2}`".format(line, n1.name, c1.name))
        
        if not hasattr(c2, n2.name):
            raise pkex.BasePyKatException("cav command `{0}`: node `{1}` is not attached to `{2}`".format(line, n2.name, c2.name))
            
        return pykat.commands.cavity(v[1], c1, n1, c2, n2)
        
        
class gauss(object):
    @staticmethod
    def parseFinesseText(text, kat):
        
        values = text.split()
        if not values[0].startswith("gauss") or (len(values) != 6 and len(values) != 8):
            raise pkex.BasePyKatException("'{0}' not a valid Finesse gauss command".format(text))        
        
        name = values[1]
        component = values[2]
        node = values[3]
        
        # setting the name of the gauss parameter is slightly convoluted
        # as we don't explicitly store gauss paramters as an object, they
        # are simply just complex numbers stored at each node. To fix this
        # the name is stored in the NodeGaussSetter object for each component
        
        if component in kat.components:
            c = kat.components[component]
            if hasattr(c, node):
                ns = getattr(c, node)
                ns.name = name
            else:
                raise pkex.BasePyKatException("Component '{0}' is not attached to node {1}".format(component, node))        
        else:
            raise pkex.BasePyKatException("Component '{0}' was not found".format(component))        
        
        if not values[0].endswith("*"):
            if len(values) == 6:
                gp = beam_param(kat.lambda0, w0=values[-2], z=values[-1])
            elif len(values) == 8:
                gpx = beam_param(kat.lambda0, w0=values[-4], z=values[-3])
                gpy = beam_param(kat.lambda0, w0=values[-2], z=values[-1])
        elif values[0].endswith("*"):
            if len(values) == 6:
                gp = beam_param(kat.lambda0, z=values[-2], zr=values[-1])
            elif len(values) == 8:
                gpx = beam_param(kat.lambda0, z=values[-4], zr=values[-3])
                gpy = beam_param(kat.lambda0, z=values[-2], zr=values[-1])
        elif values[0].endswith("**"):
            if len(values) == 6:
                gp = beam_param(kat.lambda0, w=values[-2], rc=values[-1])
            elif len(values) == 8:
                gpx = beam_param(kat.lambda0, w=values[-4], rc=values[-3])
                gpy = beam_param(kat.lambda0, w=values[-2], rc=values[-1])
        else:
            raise pkex.BasePyKatException("Unexpected ending to gauss command '{0}'".format(text))
            
        if len(values) == 6:
            kat.nodes[node].setGauss(kat.components[component], gp)
        else:
            kat.nodes[node].setGauss(kat.components[component], gpx, gpy)
            
class tf(Command):
    fQ = namedtuple('fQ', ['f', 'Q'])
    
    def __init__(self, name, poles, zeros):
        Command.__init__(self, name, False)
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
        Command.__init__(self, axis_type, True)
        
        self._axis_type = axis_type

        self.x = putter("x1")
        self.mx = putter("mx1")

        if scale == "lin":
            scale = Scale.linear
        elif scale == "log":
            scale = Scale.logarithmic
        elif isinstance(scale, six.string_types):
            # else we have a string but not a recognisable one
            raise pkex.BasePyKatException("scale argument '{0}' is not valid, must be 'lin' or 'log'".format(scale))

        if scale != Scale.linear and scale != Scale.logarithmic:
            raise pkex.BasePyKatException("scale is not Scale.linear or Scale.logarithmic")

        self.scale = scale

        if numpy.size(limits) != 2 :
            raise pkex.BasePyKatException("limits input should be a 2x1 vector of limits for the xaxis")

        self.limits = numpy.array(SIfloat(limits)).astype(float)

        if int(steps) <= 0 :
            raise pkex.BasePyKatException("steps value should be > 0")

        self.steps = int(steps)

        if isinstance(param, six.string_types):
            self.__param = param
            if comp == None:
                raise pkex.BasePyKatException("If parameter is set with a string, the comp argument must set the component name")
                
            self.__comp = comp
        elif not isinstance(param, Param) :
            raise pkex.BasePyKatException("param argument is not of type Param")
        else:
            self.__param = param
            self.__comp = param._owner()

    @property
    def param(self): return self.__param
    @param.setter
    def param(self, value):
        if not isinstance(value, Param):
            raise pkex.BasePyKatException("param argument is not of type Param")
        else:
            self.__param = value
            self.__comp = value._owner()
            
    @staticmethod
    def parseFinesseText(text):
        values = text.split()

        if values[0] != "xaxis" and values[0] != "xaxis*":
            raise pkex.BasePyKatException("'{0}' not a valid Finesse xaxis command".format(text))

        axis_type = values[0]

        values.pop(0) # remove initial value

        if len(values) != 6:
            raise pkex.BasePyKatException("xaxis Finesse code format incorrect '{0}'".format(text))

        return xaxis(values[2], [values[3], values[4]], values[1], values[5], comp=values[0], axis_type=axis_type)

    def getFinesseText(self):
        if(self._kat.noxaxis):
            return '# noaxis is true, switching xaxis off'
            
        # store either the component name of the string provided
        comp_name = self.__comp.name if hasattr(self.__comp, "name") else self.__comp
        param_name = self.__param.name if isinstance(self.__param, Param) else self.__param
        
        return '{axis_type} {0} {1} {2} {3:.16g} {4:.16g} {5}'.format(
                comp_name, param_name, self.scale,
                min(self.limits), max(self.limits), self.steps, axis_type=self._axis_type);

class x2axis(xaxis):
    def __init__(self, scale, limits, param, steps, comp=None, axis_type="x2axis"):
        xaxis.__init__(self, scale, limits, param, steps, comp=comp, axis_type=axis_type)
        self.x = putter("x2")
        self.mx = putter("mx2")

    @staticmethod
    def parseFinesseText(text):
        values = text.split()

        if values[0] != "x2axis" and values[0] != "x2axis*":
            raise pkex.BasePyKatException("'{0}' not a valid Finesse xaxis command".format(text))

        axis_type = values[0]

        values.pop(0) # remove initial value

        if len(values) != 6:
            raise pkex.BasePyKatException("xaxis Finesse code format incorrect '{0}'".format(text))

        return x2axis(values[2], [values[3], values[4]], values[1], values[5], comp=values[0],axis_type=axis_type)


class lock(Command):
    pass
