import pykat.exceptions as pkex
import numpy
import math
import copy

class gauss_param(object):
    """
    Gaussian beam complex parameter
    
    gauss_param is effectively a complex number with extra
    functionality to determine beam parameters.
    
    defaults to 1064e-9m for wavelength and refractive index 1
    usage:
        q = gauss_param(w0=w0, z=z)
        q = gauss_param(z=z, zr=zr)
        q = gauss_param(w=w, rc=rc)
        q = gauss_param(q=a) # where a is a complex number
        
        or change default wavelength and refractive index with:
        
        q = gauss_param(wavelength, nr, w0=w0, zr=zr)
    """
    
    def __init__(self, wavelength=1064e-9, nr=1, *args, **kwargs):
        self.__q = None
        self.__lambda = float(wavelength)
        self.__nr = float(nr)
        
        if len(args) == 1:
            self.__q = complex(args[0])
        
        elif len(kwargs) == 1:
            if "q" in kwargs:
                self.__q = complex(kwargs["q"])        
            else:
                raise pkex.BasePyKatException("Must specify: z and w0 or z and zr or rc and w or q, to define the beam parameter")
                
        elif len(kwargs) == 2:        
            
            if "w0" in kwargs and "z" in kwargs:
                q = float(kwargs["z"]) + 1j *float(math.pi*float(kwargs["w0"])**2/(self.__lambda/self.__nr) )
            elif "z" in kwargs and "zr" in kwargs:
                q = float(kwargs["z"]) + 1j *float(kwargs["zr"]) 
            elif "rc" in kwargs and "w" in kwargs:
                one_q = 1 / float(kwargs["rc"]) - 1j * self.__lamda / (math.pi * self.__nr * float(kwargs["w"])**2)
                q = 1/one_q
            else:
                raise pkex.BasePyKatException("Must specify: z and w0 or z and zr or rc and w or q, to define the beam parameter")
                
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
    def w(self):
        return abs(self.__q)*math.sqrt(self.__lambda / (self.__nr * math.pi * self.__q.imag))
        #return self.w0 * math.sqrt(1 + (self.__q.real/self.__q.imag)**2)
    
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
        self.__q *= complex(a)
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
        
    def __rsub__(self, a):
        return gauss_param(self.__lambda, self.__nr, complex(a) - self.__q)
    
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
        
    @property
    def real(self): return self.__q.real
    @real.setter
    def real(self, value): self.__q.real = SIfloat(value)
    
    @property
    def imag(self): return self.__q.imag
    @imag.setter
    def imag(self, value): self.__q.imag = SIfloat(value)

    # reverse beam direction 
    def reverse(self):
        self.__q = -1.0 * self.__q.real + 1j * self.__q.imag
    
class HG_gauss_beam(object):
    
    def __init__(self, qx, qy=None, n=0, m=0):
        self._qx = copy.deepcopy(qx)
        self._2pi_qrt = math.pow(2.0/math.pi, 0.25)
        
        if qy == None:
            self._q0 = copy.deepcopy(qx)
        else:
            self._q0 = copy.deepcopy(qy)
    
        self.n = n
        self.m = m
        self._calc_constants()
        
    @property
    def n(self): return self._n
    @n.setter
    def n(self,value): 
        self._n = float(value)
        self._calc_constants()

    @property
    def m(self): return self._m
    @m.setter
    def m(self,value): 
        self._m = float(value)
        self._calc_constants()
        
    def _calc_constants(self):
        self.__xpre_const = math.pow(2.0/math.pi, 0.25)
        self.__xpre_const *= math.sqrt(1/(2**self._n * math.factorial(self._n)))
        self.__xpre_const *= math.sqrt(1j*self._qx.imag / self._qx)
        self.__xpre_const *= math.pow((1j*self._qx.imag * self._qx.conjugate)/(-1j*self._qx.imag * self._qx), self._n/2.0)
        
        self.__ypre_const = math.pow(2.0/math.pi, 0.25)
        self.__ypre_const *= math.sqrt(1/(2**self._m * math.factorial(self._m)))
        self.__ypre_const *= math.sqrt(1j*self._qy.imag / self._qy)
        self.__ypre_const *= math.pow((1j*self._qy.imag * self._qy.conjugate)/(-1j*self._qy.imag * self._qy), self._m/2.0)
                
    def Unm(self, n, m, x, y):  
        # need to finish...
        return self.__xpre_const
        
        
