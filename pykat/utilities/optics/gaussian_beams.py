import pykat.exceptions as pkex
import numpy
import math

class gauss_param(object):
    """
    Gaussian beam complex parameter
    
    gauss_param is effectively a complex number with extra
    functionality to determine beam parameters.
    
    defaults to 1064e-9m for wavelength and refractive index 1
    usage:
        q = gauss_param(w0=w0, z=z)
        q = gauss_param(z=z, zr=zr)
        q = gauss_param(wz=wz, rc=rc)
        q = gauss_param(q=a) # where a is a complex number
        
        or change default wavelength and refractive index with:
        
        q = gauss_param(wavelength, nr, w0=w0, zr=zr)
    """
    
    def __init__(self, wavelength=1064e-9, nr=1, *args, **kwargs):
        self.__q = None
        self.__lambda = wavelength
        self.__nr = nr
        
        if len(args) == 1:
            self.__q = args[0]
        
        elif len(kwargs) == 1:
            if "q" in kwargs:
                self.__q = complex(kwargs["q"])        
            else:
                raise pkex.BasePyKatException("Must specify: z and w0 or z and zr or rc and wz or q, to define the beam parameter")
                
        elif len(kwargs) == 2:        
            
            if "w0" in kwargs and "z" in kwargs:
                q = float(kwargs["z"]) + 1j *float(math.pi*kwargs["w0"]**2/(self.__lambda/self.__nr) )
            elif "z" in kwargs and "zr" in kwargs:
                q = float(kwargs["z"]) + 1j *float(kwargs["zr"]) 
            elif "rc" in kwargs and "wz" in kwargs:
                one_q = 1 / kwargs["rc"] - 1j * self.__lamda / (math.pi * self.__nr * kwargs["wz"]**2)
                q = 1/one_q
            else:
                raise pkex.BasePyKatException("Must specify: z and w0 or z and zr or rc and wz or q, to define the beam parameter")
                
            self.__q = q
        else:
            raise pkex.BasePyKatException("Incorrect usage for gauss_param constructor")
    
    @property
    def wavelength(self): return self.__lambda
    
    @property
    def nr(self): return self.__nr
    
    @property
    def q(self): return self.__q
    
    @property
    def z(self): return self.__q.real
    
    @property
    def zr(self): return self.__q.imag
    
    @property
    def wz(self):
        return math.sqrt(self.__lambda /(self.__nr * math.pi) * abs(self.__q) / self.__q.imag)
    
    @property
    def w0(self):
        return math.sqrt(self.__q.imag * self.__lambda / (self.__nr * math.pi))    

    @property
    def Rc(self):
        if self.__q.real != 0:
            return abs(self.__q) / self.__q.real
        else:
            return float("inf")
    
    def conjugate(self):
        return gauss_param(self.__lambda, self.__nr, self.__q.conjugate())
    
    def __complex__(self):
        return self.__q
    
    def __str__(self):
        return str(self.__q)
    
    def __mul__(self, a):
        return gauss_param(self.__lambda, self.__nr, self.__q * complex(a))
    
    def __imul__(self, a):
        self.__q += complex(a)
        return self
        
    __rmul__ = __mul__
    
    def __add__(self, a):
        return gauss_param(self.__lambda, self.__nr, self.__q + complex(a))
    
    def __iadd__(self, a):
        self.__q += complex(a)
        return self
        
    __radd__ = __add__
    
    def __sub__(self, a):
        return gauss_param(self.__lambda, self.__nr, self.__q - complex(a))
    
    def __isub__(self, a):
        self.__q -= complex(a)
        return self
        
    __rsub__ = __sub__
    
    def __div__(self, a):
        return gauss_param(self.__lambda, self.__nr, self.__q / complex(a))
    
    def __idiv__(self, a):
        self.__q /= complex(a)
        return self
    
    def __pow__(self, q):
        return gauss_param(self.__lambda, self.__nr, self.__q**q)

    def __neg__(self, q):
        return gauss_param(self.__lambda, self.__nr, -self.__q)
        
    def __eq__(self, q):
        return complex(q) == self.__q