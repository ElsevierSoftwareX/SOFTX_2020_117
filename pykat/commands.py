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
import warnings
import pykat
import pykat.external.six as six
import pykat.exceptions as pkex

if six.PY2:
	import exceptions
    
from pykat.components import *
from pykat.structs import *
from numpy import min, max
from pykat.param import Param, putter
from collections import namedtuple
from pykat.optics.gaussian_beams import BeamParam
from pykat.freeze import canFreeze

@canFreeze
class Command(object):
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, name, unique):
        self.__dict__["____FROZEN____"] = False
        self._kat = None
        self.__unique = unique
        self.tag = None
        self.__removed = False
        self.__name = name.strip("*")
        self._putters = []
    
        
    def __deepcopy__(self, memo):
        """
        When deep copying a kat object we need to take into account
        the instance specific properties.
        """
        cls = self.__class__
        result = cls.__new__(cls)
        result._unfreeze()
        result.__dict__ = copy.deepcopy(self.__dict__, memo)
        
        for _ in result._putters:
            _._updateOwner(result)

        result._freeze()
        return result
    
    def getFinesseText(self):
        """ Base class for individual finesse optical components """
        raise NotImplementedError("This function is not implemented")

    @staticmethod
    def parseFinesseText(line, kat):
        raise NotImplementedError("This function is not implemented")

    def _on_kat_add(self, kat):
        """
        Called when this component has been added to a kat object
        """
        self._kat = kat
        
        for _ in self._putters:
            kat.registerVariable(_.name, _)

    def _on_kat_remove(self):
        self.__removed = True
        
        for i in range(len(self._putters)):
            _ = self._putters[i]
            
            self._kat.unregisterVariable(_.name)
            _.clearPuts()
        
        for i in range(len(self._putters)):  
            del self._putters[0]
            
        del self._putters[:]
        
        
    def remove(self):
        if self.__removed:
            raise pkex.BasePyKatException("{0} has already been marked as removed".format(self.name))
        else:
            self._kat.remove(self)
    
    @property
    def name(self): return self.__name
    
    @property
    def removed(self): return self.__removed





class variable(Command):
    def __init__(self, name, value):
        Command.__init__(self, name, False)
        self.__value = value
        self._freeze()
        
    def getFinesseText(self):
        return "variable {name} {value}".format(name=self.name, value=self.value)
    
    @staticmethod
    def parseFinesseText(line, kat):
        v = line.split()
        
        if len(v) != 3:
            raise pkex.BasePyKatException("'{0}' not a valid Finesse variable command".format(line))
        
        return variable(v[1], SIfloat(v[2]))
    
    @property
    def value(self): return self.__value
    @value.setter
    def value(self, Value): self.__value = SIfloat(Value)





class func(Command):
    def __init__(self, name, value):
        Command.__init__(self, name, False)
        
        self.value = value
        self.noplot = False
        self.enabled = True
        
        self.output = putter(name, self)
        self._putters.append(self.output)

        self._freeze()
        
    def getFinesseText(self):
        rtn = []

        if self.enabled:
            if self.noplot:
                rtn.append("noplot " + self.name)
        
            rtn.append("func {name} = {value}".format(name=self.name, value=str(self.value)))

        return rtn

    @staticmethod
    def parseFinesseText(line, kat):
        v = line.split(None, 3)
        v2 = line.split("=", 2)
        
        if "=" in v and len(v) == 4:
            v.remove("=")
            return func(v[1], v[2]) 
        if len(v2) == 2:
            return func(v2[0].split()[1], v2[1]) 
        else:
            raise pkex.BasePyKatException("'{0}' not a valid Finesse func command".format(line))
            

class lock(Command):
    def __init__(self, name, variable, gain, accuracy, singleLock=False):
        Command.__init__(self, name, False)
        
        self.__variable = variable
        self.__gain = gain
        self.__accuracy = accuracy
        self.singleLock = singleLock
        self.enabled = True
        self.noplot = False
        
        self.output = putter(name, self)
        self._putters.append(self.output)

        self._freeze()
        
    @staticmethod
    def parseFinesseText(line, kat):
        v = line.split()
        
        if len(v) != 5:
            raise pkex.BasePyKatException("'{0}' not a valid Finesse lock command".format(line))
            
        return lock(v[1], v[2], SIfloat(v[3]), SIfloat(v[4]), "*" in v[0])


    def getFinesseText(self):
        if self.enabled:
            cmds = "{name} {var} {gain} {accuracy}".format( name=self.name,
                                                            var=str(self.variable),
                                                            gain=str(self.gain),
                                                            accuracy=str(self.accuracy))
            
            if self.singleLock:
                return "lock* %s" % cmds
            else:
                return "lock %s" % cmds
        else:
            return None

    @property
    def variable(self): return self.__variable
    @variable.setter
    def variable(self, value): self.__variable = value

    @property
    def gain(self): return self.__gain
    @gain.setter
    def gain(self, value): self.__gain = SIfloat(value)

    @property
    def accuracy(self): return self.__accuracy
    @accuracy.setter
    def accuracy(self, value): self.__accuracy = SIfloat(value)



class cavity(Command):
    def __init__(self, name, c1, n1, c2, n2):
        Command.__init__(self, name, False)
        
        self.__c1 = c1
        self.__c2 = c2
        self.__n1 = n1
        self.__n2 = n2
        
        self.enabled = True

        self._freeze()
        
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
                gp = BeamParam(kat.lambda0, w0=values[-2], z=values[-1])
            elif len(values) == 8:
                gpx = BeamParam(kat.lambda0, w0=values[-4], z=values[-3])
                gpy = BeamParam(kat.lambda0, w0=values[-2], z=values[-1])
        elif values[0].endswith("*"):
            if len(values) == 6:
                gp = BeamParam(kat.lambda0, z=values[-2], zr=values[-1])
            elif len(values) == 8:
                gpx = BeamParam(kat.lambda0, z=values[-4], zr=values[-3])
                gpy = BeamParam(kat.lambda0, z=values[-2], zr=values[-1])
        elif values[0].endswith("**"):
            if len(values) == 6:
                gp = BeamParam(kat.lambda0, w=values[-2], rc=values[-1])
            elif len(values) == 8:
                gpx = BeamParam(kat.lambda0, w=values[-4], rc=values[-3])
                gpy = BeamParam(kat.lambda0, w=values[-2], rc=values[-1])
        else:
            raise pkex.BasePyKatException("Unexpected ending to gauss command '{0}'".format(text))
            
        if len(values) == 6:
            kat.nodes[node].setGauss(kat.components[component], gp)
        else:
            kat.nodes[node].setGauss(kat.components[component], gpx, gpy)
 
# class tf(Command):
#
#     class fQ(object):
#         def __init__(self, f, Q, tf):
#             assert(tf is not None)
#             self._tf = tf
#             self.__f = Param("f", self, None, canFsig=False, isPutable=True, isPutter=False, isTunable=True)
#             self.__Q = Param("Q", self, None, canFsig=False, isPutable=True, isPutter=False, isTunable=True)
#
#         def _register_param(self, param):
#             self._tf._params.append(param)
#
#         @property
#         def f(self): return self.__f
#         @f.setter
#         def f(self,value): self.__f.value = SIfloat(value)
#
#         @property
#         def Q(self): return self.__Q
#         @Q.setter
#         def Q(self,value): self.__Q.value = SIfloat(value)
#
#     def __init__(self, name):
#         Command.__init__(self, name, False)
#         self.zeros = []
#         self.poles = []
#         self.gain = 1
#         self.phase = 0
#         self._params = []
#
#     def addPole(self,f, Q):
#         self.poles.append(tf.fQ(SIfloat(f), SIfloat(Q), self))
#
#     def addZero(self,f, Q):
#         self.zeros.append(tf.fQ(SIfloat(f), SIfloat(Q), self))
#
#     @staticmethod
#     def parseFinesseText(text):
#         values = text.split()
#
#         if ((len(values)-4) % 3) != 0:
#             raise pkex.BasePyKatException("Transfer function Finesse code format incorrect '{0}'".format(text))
#
#         _tf = tf(values[1])
#
#         _tf.gain = SIfloat(values[2])
#         _tf.phase = SIfloat(values[3])
#
#         N = int((len(values)-4) / 3)
#
#         for i in range(1,N+1):
#             if values[i*3+1] == 'p':
#                 _tf.addPole(SIfloat(values[i*3+2]), SIfloat(values[i*3+3]))
#             elif values[i*3+1] == 'z':
#                 _tf.addZero(SIfloat(values[i*3+2]), SIfloat(values[i*3+3]))
#             else:
#                 raise pkex.BasePyKatException("Transfer function pole/zero Finesse code format incorrect '{0}'".format(text))
#
#         return _tf
#
#     def getFinesseText(self):
#         rtn = "tf {name} {gain} {phase} ".format(name=self.name,gain=self.gain,phase=self.phase)
#
#         for p in self.poles:
#             rtn += "p {f} {Q} ".format(f=p.f, Q=p.Q)
#
#         for z in self.zeros:
#             rtn += "p {f} {Q} ".format(f=z.f, Q=z.Q)
#
#         return rtn
                   
class tf(Command):
    
    class fQ(object):
        def __init__(self, f, Q):
            self.f = f
            self.Q = Q
            
    def __init__(self, name):
        Command.__init__(self, name, False)
        self.zeros = []
        self.poles = []
        self.gain = 1
        self.phase = 0
        
        self._freeze()
    
    def addPole(self,f, Q):
        self.poles.append(tf.fQ(SIfloat(f), SIfloat(Q)))
    
    def addZero(self,f, Q):
        self.zeros.append(tf.fQ(SIfloat(f), SIfloat(Q)))
    
    @staticmethod
    def parseFinesseText(text):
        values = text.split()
        
        if ((len(values)-4) % 3) != 0:
            raise pkex.BasePyKatException("Transfer function Finesse code format incorrect '{0}'".format(text))

        _tf = tf(values[1])
        
        _tf.gain = SIfloat(values[2])
        _tf.phase = SIfloat(values[3])
        
        N = int((len(values)-4) / 3)
        
        for i in range(1,N+1):
            if values[i*3+1] == 'p':
                _tf.addPole(SIfloat(values[i*3+2]), SIfloat(values[i*3+3]))
            elif values[i*3+1] == 'z':
                _tf.addZero(SIfloat(values[i*3+2]), SIfloat(values[i*3+3]))
            else:
                raise pkex.BasePyKatException("Transfer function pole/zero Finesse code format incorrect '{0}'".format(text))
    
        return _tf
        
    def getFinesseText(self):
        rtn = "tf {name} {gain} {phase} ".format(name=self.name,gain=self.gain,phase=self.phase)
        
        for p in self.poles:
            rtn += "p {f} {Q} ".format(f=p.f, Q=p.Q)
        
        for z in self.zeros:
            rtn += "p {f} {Q} ".format(f=z.f, Q=z.Q)
        
        return rtn
        
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
        
        self._set_variables()
        
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

        self._freeze()
        
    def _set_variables(self):
        self.x = putter("x1", self)
        self.mx = putter("mx1", self)

        self._putters.append(self.x)
        self._putters.append(self.mx)
        
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
                self.limits[0], self.limits[1], self.steps, axis_type=self._axis_type);

class x2axis(xaxis):
    def __init__(self, scale, limits, param, steps, comp=None, axis_type="x2axis"):
        xaxis.__init__(self, scale, limits, param, steps, comp=comp, axis_type=axis_type)

    def _set_variables(self):
        self.x = putter("x2", self)
        self.mx = putter("mx2", self)

        self._putters.append(self.x)
        self._putters.append(self.mx)
        
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

