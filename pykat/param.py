import abc
import pykat.exceptions as pkex

class putable(object):
    """
    Objects that inherit this should be able to have something `put` to it.
    Essentially this means you could write Finesse commands like
    
    put this parameter value
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, component_name, parameter_name, isPutable=True):
        self._parameter_name = parameter_name
        self._component_name = component_name
        self._putter = None
        self._isPutable  = isPutable
    
    @property
    def isPutable(self): return self._isPutable
    
    def put(self, var):
    
        if not isinstance(var, putter):
            raise pkex.BasePyKatException("var was not something that can be `put` as a value")
        
        if self._putter != None:
            self._putter.put_count -= 1
        
        self._putter = var
        self._putter.put_count += 1
        
    def _getPutFinesseText(self):
        rtn = []
        # if something is being put to this 
        if self._putter != None:
            rtn.append("put {comp} {param} ${value}".format(comp=self._component_name, param=self._parameter_name, value=self._putter.put_name()))
        
        return rtn
            
class putter(object):
    """
    If an object can be put to something that is putable it should inherit this
    object.
    """
    
    def __init__(self, put_name, isPutter=True):
        self._put_name = put_name
        self.put_count = 0
        self._isPutter = isPutter
    
    @property
    def isPutter(self): return self._isPutter
    
    def put_name(self): return self._put_name
    
        
class Param(putable, putter):

    def __init__(self, name, owner, value, isPutable=True, isPutter=True, isTunable=True, var_name=None):
        self._name = name
        self._owner = owner
        self._value = value
        self._isPutter = isPutter
        self._isTunable = isTunable
        self._owner._register_param(self)
        
        if isPutter:
            if var_name == None:
                var_name = "var_{0}_{1}".format(owner.name, name)
                
            putter.__init__(self, var_name, isPutter)
            
        if isPutable:
            putable.__init__(self, owner.name, name, isPutable)
        
    @property
    def name(self): return self._name
    
    @property
    def isTuneable(self): return self._isTunable
    
    @property
    def value(self): return self._value
    @value.setter
    def value(self, value):
        self._value = value
    
    def __str__(self): return str(self.value)
    def __float__(self): return self.value
        
    def getFinesseText(self):
        rtn = []
        
        if self.isPutable: rtn.extend(self._getPutFinesseText())
        
        # if this parameter is being put somewhere then we need to
        # set it as a variable
        if self.isPutter and self.put_count > 0:
            rtn.append("set {put_name} {comp} {param}".format(put_name=self.put_name(), comp=self._owner.name, param=self.name))
        
        return rtn
        
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
        
    __rsub__ = __sub__
    
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
        rtn = []
        
        if self.value != 0:
            rtn.append("attr {0} {1} {2}".format(self._owner.name, self.name, self.value))
            
        rtn.extend(super(AttrParam, self).getFinesseText())
        
        return rtn

    
