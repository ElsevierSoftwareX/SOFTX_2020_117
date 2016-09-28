from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import abc
import pykat.exceptions as pkex
import weakref
    
class putable(object):
    """
    Objects that inherit this should be able to have something `put` to it.
    Essentially this means you could write Finesse commands like
    
    param.put(kat.xaxis.x)
    
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, component_name, parameter_name, isPutable=True):
        self._parameter_name = parameter_name
        self._component_name = component_name
        self._putter = None
        self._alt = False
        self._isPutable  = isPutable
    
    @property
    def isPutable(self): return self._isPutable
    
    def put(self, var, alt=False):
        if not self._isPutable:
            raise pkex.BasePyKatException("Can't put to this object")
            
        if var is not None and not isinstance(var, putter):
            raise pkex.BasePyKatException("`%s` was not something that can be `put` to a parameter" % str(var))
        
        # Remove existing puts
        if self._putter is not None:
            self._putter.unregister(self)
        
        self._putter = var
            
        self._alt = alt
        
        if var is not None:
            self._putter.register(self)
        
    def _getPutFinesseText(self):
        rtn = []

        if self._isPutable and self._putter is not None:
            putter_enabled = True
                
            if hasattr(self._putter.owner, 'enabled'):
                putter_enabled = self._putter.owner.enabled
                                
            if putter_enabled:
                if self._alt:
                    alt = '*'
                else:
                    alt = ''
    
                # if something is being put to this 
                rtn.append("put{alt} {comp} {param} ${value}".format(alt=alt, comp=self._component_name, param=self._parameter_name, value=self._putter.put_name()))
        
        return rtn
        
        
class putter(object):
    """
    If an object can be put to something that is putable it should inherit this
    object.
    """
    
    def __init__(self, put_name, owner, isPutter=True):
        self._put_name = put_name
        self.put_count = 0
        self._isPutter = isPutter
        self.putees = [] # list of params that this puts to
        
        assert(owner is not None)
        self.__owner = weakref.ref(owner)
    
    def _updateOwner(self, newOwner):
        del self.__owner
        self.__owner = weakref.ref(newOwner)
        
    def clearPuts(self):
        import copy
        for _ in copy.copy(self.putees):
            _.put(None)
    
    def register(self, toput):
        if not self._isPutter:
            raise pkex.BasePyKatException("This object can't put")
            
        self.put_count += 1
        self.putees.append(toput)
    
    def unregister(self, item):
        if not self._isPutter:
            raise pkex.BasePyKatException("This object can't put")
            
        self.put_count -= 1
        self.putees.remove(item)
        
    @property
    def owner(self): return self.__owner()
    
    @property
    def name(self): return self._put_name
    
    @property
    def putCount(self): return self.put_count
    
    @property
    def isPutter(self): return self._isPutter
    
    def put_name(self): return self._put_name
    
        
class Param(putable, putter):

    def __init__(self, name, owner, value, canFsig=False, fsig_name=None, isPutable=True, isPutter=True, isTunable=True, var_name=None, register=True):
        self._name = name
        self._registered = register
        self._owner = weakref.ref(owner)
        self._value = value
        self._isPutter = isPutter
        self._isTunable = isTunable
        self._canFsig = False
        
        if self._registered:
            self._owner()._register_param(self)
        
        if canFsig:
            self._canFsig = True
            
            if fsig_name is None:
                raise pkex.BasePyKatException("If parameter is a possible fsig target the fsig_name argument must be set")
                
            self.__fsig_name = fsig_name
        
        if isPutter:
            if var_name is None:
                var_name = "var_{0}_{1}".format(owner.name, name)
                
        putter.__init__(self, var_name, owner, isPutter)
            
        putable.__init__(self, owner.name, name, isPutable)
        
    @property
    def canFsig(self): return self._canFsig
    
    @property
    def owner(self): return self._owner()
    
    @property
    def fsig_name(self): return self.__fsig_name
    
    @property
    def fsigName(self): return self.__fsig_name
    
    @property
    def name(self): return self._name
    
    @property
    def isTuneable(self): return self._isTunable
    
    @property
    def value(self):
        if self._owner().removed:
            raise pkex.BasePyKatException("{0} has been removed from the simulation".format(self._owner().name))
        else:
            return self._value
    
    @value.setter
    def value(self, value):
        if self._owner().removed:
            raise pkex.BasePyKatException("{0} has been removed from the simulation".format(self._owner().name))
        else:
            self._value = value
    
    def __str__(self):
        if self._owner().removed:
            raise pkex.BasePyKatException("{0} has been removed from the simulation".format(self._owner().name))
        elif type(self.value) == float:
            return repr(self.value)
        else:
            return str(self.value)
            
    def __float__(self):
        if self._owner().removed:
            raise pkex.BasePyKatException("{0} has been removed from the simulation".format(self._owner().name))
        else:
            return float(self.value)
        
    def getFinesseText(self):
        if self._owner().removed:
            raise pkex.BasePyKatException("{0} has been removed from the simulation".format(self._owner().name))
            
        rtn = []
        
        if self.isPutable: rtn.extend(self._getPutFinesseText())
        
        # if this parameter is being put somewhere then we need to
        # set it as a variable
        if self.isPutter and self.put_count > 0:
            rtn.append("set {put_name} {comp} {param}".format(put_name=self.put_name(), comp=self._owner().name, param=self.name))
        
        return rtn
        
    def _updateOwner(self, newOwner):
        """
        This updates the internal weak reference to link a parameter to who owns it.
        Should only be called by the __deepcopy__ component method to ensure things
        are kept up to date.
        """
        del self._owner
        self._owner = weakref.ref(newOwner)
        
    def _onOwnerRemoved(self):
        #if this param can be put somewhere we need to check if it is
        if self.isPutable:
            for a in self.putees:
                print ("Removing put from {0} {1} to {2} {3}".format(self.owner.name, self.name, a.owner.name, a.name))
                a._putter = None
                self.put_count -= 1
                
            # delete any references left over
            del self.putees[:]
        
        # check if we have anything being put to us
        if self.isPutter:
            if self._putter != None:
                print ("Removing put from {0} {1} to {2} {3}".format(self._putter.owner.name, self._putter.name, self.owner.name, self.name))
                self._putter.put_count -= 1
                self._putter.putees.remove(self)
                self._putter = None
       
       
    def __mul__(self, a):
        return self.value * a
    
    def __imul__(self, a):
        return self.value * (a)
        
    __rmul__ = __mul__
    
    def __add__(self, a):
        return self.value + (a)
    
    def __iadd__(self, a):
        return self.value + (a)
        
    __radd__ = __add__
    
    def __sub__(self, a):
        return self.value - (a)
    
    def __isub__(self, a):
        return self.value - (a)
        
    def __rsub__(self, a):
        return (a) - self.value
    
    def __div__(self, a):
        return self.value / (a)
    
    def __idiv__(self, a):
        return self.value / complex(a)
        
    def __pow__(self, q):
        return  self.value**q

    def __neg__(self):
        return -self.value
        
    def __eq__(self, q):
        return (q) == self.value
    def __ne__(self, q):
        return (q) != self.value
    def __lt__(self, q):
        return (q) > self.value
    def __gt__(self, q):
        return (q) < self.value        
        
class AttrParam(Param):
    """
    Certain parameters of a component are set using the Finesse `attr` command.
    
    This inherits directly from a Param object so can be set whether this attribute
    is putable or a putter.
    
    If the value pf the parameter is not 0 the attr command will be printed.
    """
    def getFinesseText(self):
        if self._owner().removed:
            raise pkex.BasePyKatException("{0} has been removed from the simulation".format(self._owner().name))

        rtn = []
        
        if self.value != None:
            rtn.append("attr {0} {1} {2}".format(self._owner().name, self.name, self.value))
            
        rtn.extend(super(AttrParam, self).getFinesseText())
        
        return rtn

    
