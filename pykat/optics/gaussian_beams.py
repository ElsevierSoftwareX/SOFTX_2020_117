from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pykat.exceptions as pkex
import numpy as np
import math
import copy
import warnings
import cmath
from math import factorial
from scipy.special import hermite
from pykat.math.jacobi import jacobi
from pykat.SIfloat import SIfloat

        
class BeamParam(object):
    """
    Gaussian beam complex parameter.
    
    BeamParam is effectively a complex number with extra
    functionality to determine beam parameters.
    
    Defaults to 1064e-9m for wavelength and refractive index 1
    usage:
        q = BeamParam(w0=w0, z=z)
        q = BeamParam(z=z, zr=zr)
        q = BeamParam(w=w, rc=rc)
        q = BeamParam(q=a) # where a is a complex number
        
        or change default wavelength and refractive index with:
        
        q = BeamParam(wavelength, nr, w0=w0, zr=zr)
    """
    
    def __init__(self, wavelength=1064e-9, nr=1, *args, **kwargs):
        if self.__class__ != BeamParam:
            warnings.warn("Name changed. Use BeamParam instead of gauss_param or beam_param.")
            
        self.__q = None
        self.__lambda = SIfloat(wavelength)
        self.__nr = SIfloat(nr)
        
        if len(args) == 1:
            self.__q = complex(args[0])
        
        elif len(kwargs) == 1:
            if "q" in kwargs:
                self.__q = complex(kwargs["q"])        
            else:
                raise pkex.BasePyKatException("Must specify: z and w0 or z and zr or rc and w or q, to define the beam parameter")
                
        elif len(kwargs) == 2:        
            
            if "w0" in kwargs and "z" in kwargs:
                q = SIfloat(kwargs["z"]) + 1j * math.pi*SIfloat(kwargs["w0"])**2/(self.__lambda/self.__nr)
            elif "z" in kwargs and "zr" in kwargs:
                q = SIfloat(kwargs["z"]) + 1j * SIfloat(kwargs["zr"]) 
            elif "rc" in kwargs and "w" in kwargs:
                one_q = 1 / SIfloat(kwargs["rc"]) - 1j * SIfloat(wavelength) / (math.pi * SIfloat(nr) * SIfloat(kwargs["w"])**2)
                q = 1/one_q
            else:
                raise pkex.BasePyKatException("Must specify: z and w0 or z and zr or rc and w or q, to define the beam parameter")
                
            self.__q = q
        else:
            raise pkex.BasePyKatException("Incorrect usage for gauss_param constructor")
    
    @property
    def wavelength(self): return self.__lambda
    @wavelength.setter
    def wavelength(self,value): self.__lambda = SIfloat(value)
    
    @property
    def nr(self): return self.__nr
    
    @property
    def q(self): return self.__q
    
    @property
    def z(self): return self.__q.real
    
    @property
    def zr(self): return self.__q.imag
    
    @property
    def w(self, z=None):
        return np.abs(self.__q)* np.sqrt(self.__lambda / (self.__nr * math.pi * self.__q.imag))
    
    def beamsize(self, z=None, wavelength=None, nr=None, w0=None):

        if z is None:
            z = self.z
        else:
            z = np.array(z)
                
        if wavelength is None:
            wavelength = self.wavelength
        else:
            wavelength = np.array(wavelength)
            
        if nr is None:
            nr = self.nr
        else:
            nr = np.array(nr)
            
        if w0 is None:
            w0 = self.w0
        else:
            w0 = np.array(w0)
        
        q = z + 1j * math.pi * w0 **2 / wavelength
        
        return np.abs(q)*np.sqrt(wavelength / (nr * math.pi * q.imag))
    
    def gouy(self, z=None, wavelength=None, nr=None, w0=None):
        if z is None:
            z = self.z
        else:
            z = np.array(z)
                
        if wavelength is None:
            wavelength = self.wavelength
        else:
            wavelength = np.array(wavelength)
            
        if nr is None:
            nr = self.nr
        else:
            nr = np.array(nr)
            
        if w0 is None:
            w0 = self.w0
        else:
            w0 = np.array(w0)
        
        q = z + 1j * math.pi * w0 **2 / wavelength
        
        return np.arctan2(q.real, q.imag)
        
    @property
    def w0(self):
        return np.sqrt(self.__q.imag * self.__lambda / (self.__nr * math.pi))    

    @property
    def Rc(self):
        def __rc(z, zr):
            if z != 0:
                return z * (1 + (zr/z)**2)
            else:
                return float("inf")
                
        v = np.vectorize(__rc)
        
        return v(self.z, self.zr)
    
    def curvature(self, z=None, wavelength=None, nr=None, w0=None):
        if z is None:
            z = self.z
        else:
            z = np.array(z)
                
        if wavelength is None:
            wavelength = self.wavelength
        else:
            wavelength = np.array(wavelength)
            
        if nr is None:
            nr = self.nr
        else:
            nr = np.array(nr)
            
        if w0 is None:
            w0 = self.w0
        else:
            w0 = np.array(w0)
        
        q = z + 1j * math.pi * w0 **2 / wavelength
        
        return q.real * (1+ (q.imag/q.real)**2)
        
    @staticmethod
    def overlap(q1, q2):
        """
        Computes the projection from one beam parameter to another to give a measure of the
        overlap between the two beam parameters.
        
        This function was provided by Paul Fulda and Antonio Perreca, which came originally
        from Chris Mueller.
        
        Added on 20/4/2015
        """
        return abs(4*q1.imag * q2.imag)/abs(q1.conjugate()-q2)**2
        
    def conjugate(self):
        return BeamParam(self.__lambda, self.__nr, self.__q.conjugate())
    
    def __abs__(self):
        return abs(complex(self.__q))
        
    def __complex__(self):
        return self.__q
    
    def __str__(self):
        return str(self.__q)
    
    def __mul__(self, a):
        return BeamParam(self.__lambda, self.__nr, self.__q * complex(a))
    
    def __imul__(self, a):
        self.__q *= complex(a)
        return self
        
    __rmul__ = __mul__
    
    def __add__(self, a):
        return BeamParam(self.__lambda, self.__nr, self.__q + complex(a))
    
    def __iadd__(self, a):
        self.__q += complex(a)
        return self
        
    __radd__ = __add__
    
    def __sub__(self, a):
        return BeamParam(self.__lambda, self.__nr, self.__q - complex(a))
    
    def __isub__(self, a):
        self.__q -= complex(a)
        return self
        
    def __rsub__(self, a):
        return BeamParam(self.__lambda, self.__nr, complex(a) - self.__q)
    
    def __div__(self, a):
        return BeamParam(self.__lambda, self.__nr, self.__q / complex(a))
    
    def __idiv__(self, a):
        self.__q /= complex(a)
        return self
    
    def __pow__(self, q):
        return BeamParam(self.__lambda, self.__nr, self.__q**q)

    def __neg__(self, q):
        return BeamParam(self.__lambda, self.__nr, -self.__q)
        
    def __eq__(self, q):
        if q is None:
            return False
            
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
        
        
class HG_mode(object):
    """ Hermite-Gauss mode profile. Example usage:
    import pykat.optics.gaussian_beams as gb
    qx=gb.BeamParam(w0=1e-3,z=0)
    beam=gb.HG_mode(qx,n=2,m=0)
    beam.plot()
    """    
    def __init__(self, qx, qy=None, n=0, m=0):
        self._qx = copy.deepcopy(qx)
        self._2pi_qrt = math.pow(2.0/math.pi, 0.25)
        
        if qy is None:
            self._qy = copy.deepcopy(qx)
        else:
            self._qy = copy.deepcopy(qy)
    
        self._n = int(n)
        self._m = int(m)
        self._hn = hermite(self._n)
        self._hm = hermite(self._m)
        self._calc_constants()
        
    @property
    def n(self): return self._n
    @n.setter
    def n(self,value): 
        self._n = int(value)
        self._calc_constants()
        self._hn = hermite(self._n)

    @property
    def m(self): return self._m
    @m.setter
    def m(self,value): 
        self._m = int(value)
        self._calc_constants()
        self._hm = hermite(self._m)
            
    @property
    def q(self):
        if self._qx.q == self._qy.q:
            return self._qx.q
        else:
            return (self._qx.q, self._qy.q)
    @q.setter
    def q(self, value):
        if value.__class__ == BeamParam:
            self._qx = copy.deepcopy(value)
            self._qy = copy.deepcopy(value)
        else:
            self._qx = BeamParam(q=complex(value))
            self._qy = BeamParam(q=complex(value))
    
    @property
    def qx(self):
        return self._qx.q
        
    @qx.setter
    def qx(self, value):
        if value.__class__ == BeamParam:
            self._qx = copy.deepcopy(value)
        else:
            self._qx = BeamParam(q=complex(value))
    
    @property
    def qy(self):
        return self._qy.q
        
    @qy.setter
    def qy(self, value):
        if value.__class__ == BeamParam:
            self._qy = copy.deepcopy(value)
        else:
            self._qy = BeamParam(q=complex(value))
    
    @property
    def constant_x(self):
        return self.__xpre_const
        
    @property
    def constant_y(self):
        return self.__ypre_const
        
    def _calc_constants(self):
        self.__xpre_const = math.pow(2.0/math.pi, 0.25)
        self.__xpre_const *= np.sqrt(1.0/(self._qx.w0 * 2**(self._n) * np.math.factorial(self._n)))
        self.__xpre_const *= np.sqrt(1j*self._qx.imag / self._qx.q)
        self.__xpre_const *= ((1j*self._qx.imag * self._qx.q.conjugate())/(-1j*self._qx.imag * self._qx.q)) ** ( self._n/2.0)
        
        self.__ypre_const = math.pow(2.0/math.pi, 0.25)
        self.__ypre_const *= np.sqrt(1.0/(self._qy.w0 * 2**(self._m) * np.math.factorial(self._m)))
        self.__ypre_const *= np.sqrt(1j*self._qy.imag / self._qy.q)
        self.__ypre_const *= ((1j*self._qy.imag * self._qy.q.conjugate())/(-1j*self._qy.imag * self._qy.q)) **(self._m/2.0)
    
        self.__sqrt2_wxz = math.sqrt(2) / self._qx.w
        self.__sqrt2_wyz = math.sqrt(2) / self._qy.w
        
        self.__kx =  2*math.pi / self._qx.wavelength
        self.__ky =  2*math.pi / self._qy.wavelength
        
        self.__invqx = 1/ self._qx.q
        self.__invqy = 1/ self._qy.q
    
    def Un(self, x):
        return self.__xpre_const * self._hn(self.__sqrt2_wxz * x) * np.exp(-0.5j * self.__kx * x*x * self.__invqx)
    
    def Um(self, y):
        return self.__ypre_const * self._hm(self.__sqrt2_wyz * y) * np.exp(-0.5j * self.__ky * y*y * self.__invqy)
        
    def Unm(self, x, y):
        _un = self.Un(x)  
        _um = self.Um(y)
        return np.outer(_un, _um)
        
    def plot(self, ndx=100, ndy=100, xscale=4, yscale=4):
        """ Make a simple plot the HG_mode """
        import pykat.plotting 
        import matplotlib.pyplot as plt
        
        xrange = xscale * np.linspace(-self._qx.w, self._qx.w, ndx)
        yrange = yscale * np.linspace(-self._qy.w, self._qy.w, ndy)

        dx = xrange[1]-xrange[0]
        dy = yrange[1]-yrange[0]

        data = self.Unm(xrange,yrange)

        fig = pykat.plotting.figure()
        axes = plt.imshow(np.abs(data.T), aspect=dx/dy, extent=[min(xrange),max(xrange),min(yrange),max(yrange)])
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        cbar = fig.colorbar(axes)
        plt.show()
        

def HG2LG(n,m):
    """A function for Matlab which returns the coefficients and mode indices of
    the LG modes required to create a particular HG mode.
    Usage: coefficients,ps,ls = HG2LG(n,m)
    
    n,m:          Indces of the HG mode to re-create.
    coeffcients:  Complex coefficients for each order=n+m LG mode required to
                  re-create HG_n,m.
    ps,ls:        LG mode indices corresponding to coefficients.
    """
    # Mode order
    N = n+m;
    
    # Create empty vectors for LG coefficients/ indices
    coefficients = np.linspace(0,0,N+1,dtype=np.complex_)
    ps = np.linspace(0,0,N+1)
    ls = np.linspace(0,0,N+1)
    
    # Calculate coefficients
    for j in np.arange(0,N+1):
        
        # Indices for coefficients
        l = 2*j-N
        p = int((N-np.abs(l))/2)
        
        ps[j] = p
        ls[j] = l
        
        signl = np.sign(l)
        if (l==0):
            signl = 1.0

        # Coefficient
        c = (signl*1j)**m * np.sqrt(factorial(N-m)*factorial(m)/(2**N * factorial(np.abs(l)+p)*factorial(p)))
        c = c * (-1.0)**p * (-2.0)**m * jacobi(m,np.abs(l)+p-m,p-m,0.0)

        coefficients[j] = c
        
    return coefficients, ps, ls 
        

    
def LG2HG(p,l):
    """ Function to compute the amplitude coefficients
    of Hermite-Gauss modes whose sum yields a Laguerre Gauss mode
    of order n,m.
    Usage: coefficients, ns, ms = LG2HG(p,l)
    p:    Radial LG index
    l:    Azimuthal LG index
    The LG mode is written as u_pl with 0<=|l|<=p.
    The output is a series of coefficients for u_nm modes,
    given as complex numbers and respective n,m indices
    coefficients (complex array): field amplitude for mode u_nm
    ns (int array): n-index of u_nm
    ms (int array): m-index of u_nm

    
    The formula is adpated from M.W. Beijersbergen et al 'Astigmatic
    laser mode converters and transfer of orbital angular momentum',
    Opt. Comm. 96 123-132 (1993)
    We adapted coefficients to be compatible with our
    definition of an LG mode, it differs from
    Beijersbergen by a (-1)^p factor and has exp(il\phi) rather
    than exp(-il\phi).  Also adapted for allowing -l.
    Andreas Freise, Charlotte Bond    25.03.2007"""

    # Mode order
    N=2*p+np.abs(l)

    # Indices for coefficients
    n = np.abs(l)+p
    m = p

    # create empty vectors
    coefficients = np.linspace(0,0,N+1,dtype=np.complex_)
    ns = np.linspace(0,0,N+1)
    ms = np.linspace(0,0,N+1)

    # l positive or negative
    signl = np.sign(l)
    if (l==0):
        signl = 1.0
    
    # Beijersbergen coefficients
    for j in np.arange(0,N+1):
        ns[j]=N-j
        ms[j]=j

        c=(-signl*1j)**j * math.sqrt(factorial(N-j)*factorial(j)/(2**N * factorial(n)*factorial(m)))
        coefficients[j] = c * (-1.0)**p * (-2)**j * jacobi(j,n-j,m-j,0.0)
    
    return coefficients, ns, ms


# These classes are here as legacy classes, BeamParam should throw a warning if they are used instead.


class gauss_param(BeamParam):
    pass

class beam_param(BeamParam):
    pass