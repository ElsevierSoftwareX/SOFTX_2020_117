import warnings

def canFreeze(cls):
    
    def _freeze(self): self.__dict__["____FROZEN____"] = True
    def _unfreeze(self): self.__dict__["____FROZEN____"] = False
    
    def frozensetattr(self, name, value):
        if "____FROZEN____" in self.__dict__ and self.__dict__["____FROZEN____"] and not hasattr(self, name):
            if hasattr(self, "name"):
                n = self.name
            elif hasattr(self, "__name"):
                n = self.__name
            else:
                n = self.__class__.__name__
                
            warnings.warn("'%s' does not have attribute called '%s'" % (n, name), stacklevel=2)
            
        super(cls, self).__setattr__(name, value)

    cls.__setattr__ = frozensetattr
    cls._freeze = _freeze
    cls._unfreeze = _unfreeze
    
    return cls