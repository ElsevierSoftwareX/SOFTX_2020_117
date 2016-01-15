import pykat.exceptions as pkex
import numpy as np
import math
import copy
import warnings
import cmath
from scipy.special import hermite
from pykat.SIfloat import SIfloat

def fit_beam_size(z_data, w_data, w0guess, z0guess, lam=1064e-9, show_plot=True,
                  plotpts=100, title='Beam scan fit', w2zero=True, filename=None):
    """
    Function to fit a beam size (w0 and z) to a 2D intensity array. Plotting can
    be switched on or off.

    Returns best fit of (w0, z)
    
    Contributed by Paul Fulda on 15/01/16
    """
    
    import scipy.optimize
    
    def zw0z02w(z, w0, z0, lam):
        z_R = pl.pi*w0**2/lam
        w = w0*pl.sqrt(1+((z-z0)/z_R)**2)
        return w
    
    popt,pcov = scipy.optimize.curve_fit(lambda z_data, w0, z0: zw0z02w(z_data, w0, z0, lam), z_data, w_data, p0=[w0guess,z0guess])
    w0out=popt[0]
    z0out=popt[1]
    
    z_fit = pl.linspace(min(z_data),max(z_data),plotpts)
    w_fit = zw0z02w(z_fit, w0out, z0out, lam)


    if show_plot or filename is not None:
        import pylab as pl
        
        um=1e6
        
        pl.figure()
        pl.plot(z_data,w_data*um,'bo', label = 'Data')
        pl.plot(z_fit,w_fit*um,'b', label = 'Fit')
        pl.tight_layout
        pl.grid()
        pl.xlabel('Position [m]')
        pl.ylabel('Beam size [$\mu$m]')
        pl.xlim((min(z_data),max(z_data)))
        
        if w2zero:
            pl.ylim((0,max(w_data)*um))
            
        pl.legend(loc=0)
        pl.title(title)
        
        if filename is not None:
            pl.savefig(filename)
            
    return w0out, z0out      
    
    

class gauss_param(object):
    """
    Use beam_param instead, will be the future name of this object.
    
    Gaussian beam complex parameter
    
    beam_param is effectively a complex number with extra
    functionality to determine beam parameters.
    
    Defaults to 1064e-9m for wavelength and refractive index 1
    usage:
        q = gauss_param(w0=w0, z=z)
        q = gauss_param(z=z, zr=zr)
        q = gauss_param(w=w, rc=rc)
        q = gauss_param(q=a) # where a is a complex number
        
        or change default wavelength and refractive index with:
        
        q = gauss_param(wavelength, nr, w0=w0, zr=zr)
    """
    
    def __init__(self, wavelength=1064e-9, nr=1, *args, **kwargs):
        if self.__class__ != beam_param:
            warnings.warn("Name changed. Use beam_param instead of gauss_param.")
            
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
    def w(self):
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
        return beam_param(self.__lambda, self.__nr, self.__q.conjugate())
    
    def __abs__(self):
        return abs(complex(self.__q))
        
    def __complex__(self):
        return self.__q
    
    def __str__(self):
        return str(self.__q)
    
    def __mul__(self, a):
        return beam_param(self.__lambda, self.__nr, self.__q * complex(a))
    
    def __imul__(self, a):
        self.__q *= complex(a)
        return self
        
    __rmul__ = __mul__
    
    def __add__(self, a):
        return beam_param(self.__lambda, self.__nr, self.__q + complex(a))
    
    def __iadd__(self, a):
        self.__q += complex(a)
        return self
        
    __radd__ = __add__
    
    def __sub__(self, a):
        return beam_param(self.__lambda, self.__nr, self.__q - complex(a))
    
    def __isub__(self, a):
        self.__q -= complex(a)
        return self
        
    def __rsub__(self, a):
        return beam_param(self.__lambda, self.__nr, complex(a) - self.__q)
    
    def __div__(self, a):
        return beam_param(self.__lambda, self.__nr, self.__q / complex(a))
    
    def __idiv__(self, a):
        self.__q /= complex(a)
        return self
    
    def __pow__(self, q):
        return beam_param(self.__lambda, self.__nr, self.__q**q)

    def __neg__(self, q):
        return beam_param(self.__lambda, self.__nr, -self.__q)
        
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


class beam_param(gauss_param):
    pass
    
class HG_beam(object):
    
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
        if value.__class__ == beam_param:
            self._qx = copy.deepcopy(value)
            self._qy = copy.deepcopy(value)
        else:
            self._qx = beam_param(q=complex(value))
            self._qy = beam_param(q=complex(value))
    
    @property
    def qx(self):
        return self._qx.q
        
    @qx.setter
    def qx(self, value):
        if value.__class__ == beam_param:
            self._qx = copy.deepcopy(value)
        else:
            self._qx = beam_param(q=complex(value))
    
    @property
    def qy(self):
        return self._qy.q
        
    @qy.setter
    def qy(self, value):
        if value.__class__ == beam_param:
            self._qy = copy.deepcopy(value)
        else:
            self._qy = beam_param(q=complex(value))
    
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
        import pylab
        
        xrange = xscale * np.linspace(-self._qx.w, self._qx.w, ndx)
        yrange = yscale * np.linspace(-self._qy.w, self._qy.w, ndy)

        dx = xrange[1]-xrange[0]
        dy = yrange[1]-yrange[0]

        data = self.Unm(xrange,yrange)

        fig = pylab.figure()
        axes = pylab.imshow(np.abs(data), aspect=dx/dy, extent=[min(xrange),max(xrange),min(yrange),max(yrange)])
        pylab.xlabel('x [m]')
        pylab.ylabel('y [m]')
        cbar = fig.colorbar(axes)
        pylab.show()
        
        
