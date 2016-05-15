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
from scipy.integrate import quad
from scipy.special import hermite
from pykat.SIfloat import SIfloat
from pykat.optics.gaussian_beams import HG2LG
from pykat.math.jacobi import jacobi
from pykat.math.laguerre import laguerre


def finesse_split_photodiode(maxtem, xy='x'):
    """Prints beat coefficients for HG modes on a split photo detector
    in the format for insertion into the kat.ini file of Finesse.
    Example use:
    finesse_split_photodiode(5, xy='x')
    The default kat.ini numbers have been generated with:
    finesse_split_photodiode(40, xy='x')
    finesse_split_photodiode(40, xy='y')
    """
    assert (xy=="x" or xy=="y"), "xy function parameter must be 'x' or 'y'" 
    assert (isinstance(maxtem, int) and maxtem >=0), "maxtem must be integer >=0" 

    if (xy=="x"):
        print('PDTYPE x-split')
    else:
        print('PDTYPE y-split')

    # running through even modes
    for n1 in np.arange(0,maxtem+1,2):
        # running through odd modes
        for n2 in np.arange(1,maxtem+1,2) :
            c_n1n2= HG_split_diode_coefficient_numerical(n1,n2)
            if (xy=="x"):       
                print("{:2d} x {:2d} x {:.15g}".format(n1,n2,c_n1n2))
            else:
                print("x {:2d} x {:2d} {:.15g}".format(n1,n2,c_n1n2))
    print('END')

                
def HG_split_diode_coefficient_analytical(n1,n2):
    """ Using an analytical equation to compute beat coefficients
    betwween mode n1 and n2 on a split photo detector. Uses arbitrary
    precision (from gmpy) because otherwise errors get very large
    for n1,n2>30. This is for comparison with the numerical method
    HG_split_diode_coefficient_numerical only. """
    import gmpy
    
    temp=gmpy.mpq(0.0)
    for l in np.arange(0,n1/2+1):
        for k in np.arange(0,(n2-1)/2+1):
            temp +=  gmpy.mpq(pow((-1.0/4.0),l+k)*gmpy.mpz(factorial((n2+n1-1)/2-l-k))/gmpy.mpz((factorial(l)*factorial(k)*factorial(n1-2*l)*factorial(n2-2*k))))

    c_n1n2=gmpy.mpq(temp*gmpy.mpq(math.sqrt(2.0**(n1+n2)*gmpy.mpz(factorial(n1)*factorial(n2))/np.pi)))
    return float(gmpy.mpf(c_n1n2))
    
    
def HG_split_diode_coefficient_numerical(n1,n2):
    """ Compute beam coefficient between mode n1 and n2 on a split
    photo detector using numerical integration, tested up to n1,n2 = 40.
    This is primarily used by finesse_split_photodiode()"""
    A = 2.0 * np.sqrt(2.0/np.pi) * np.sqrt(1.0 / (2.0**(n1+n2) * factorial(n1) * factorial(n2)))
    f = lambda x: hermite(n1)(np.sqrt(2.0)*x) * math.exp(-2.0*x*x) * hermite(n2)(np.sqrt(2.0)*x)    
    c_n1n2= quad(f, 0.0, np.Inf, epsabs=1e-10, epsrel=1e-10, limit=200)
    return A*c_n1n2[0]
    

def LG_bullseye_coefficients_numerical(l1, p1, l2, p2, w, r):
    """ Function to compute a beat coefficient for two LG modes on a
    bullseye photo detector, used by finesse_bullseye_photodiode().
    l1,p1 mode indices of first LG mode
    l2,p2 mode indices of second LG mode
    w: beam radius on diode [m]
    r: radius of inner disk element [m]
    w = beam radisu [m]
    r = radius of inner element
    """
    if (l1 != l2):
        return 0.0
        
    l = np.abs(l1)
    # computer pre-factor of integral
    A = math.sqrt(factorial(p1)*factorial(p2)/(factorial(l+p1)*factorial(l+p2)))

    # definng integrand
    f = lambda x: x**l * laguerre(p1,l,x) * laguerre(p2,l,x) * math.exp(-x)
    # defining integration limit
    int_limit = 2.0 * r**2 / w**2
    # perform numerical integrations
    val1, res= quad(f, 0.0, int_limit, epsabs=1e-10, epsrel=1e-10, limit=500)
    val2, res= quad(f, int_limit, np.Inf, epsabs=1e-10, epsrel=1e-10, limit=500)
    return A*(val2-val1)

def HG_bullseye_coefficients_numerical(n1, m1, n2, m2, w, r):
    """ Function to compute a beat coefficient for two HG modes on a
    bullseye photo detector, used by finesse_bullseye_photodiode().
    n1,m1 mode indices of first HG mode
    n2,m2 mode indices of second HG mode
    w: beam radius on diode [m]
    r: radius of inner disk element [m]

    """
    mparity=np.mod(m1+m2,2)
    nparity=np.mod(n1+n2,2)
                
    k = 0.0 + 0.0*1j
    if (mparity==0 and nparity==0):            
        a,p1,l1 = HG2LG(n1,m1)
        b,p2,l2 = HG2LG(n2,m2)
        for idx1, aa in enumerate(l1):
            for idx2, bb in enumerate(l2):
                if l1[idx1]==l2[idx2]:
                    c = LG_bullseye_coefficients_numerical(l1[idx1],p1[idx1], l2[idx2], p2[idx2], w, r)
                    k = k + a[idx1] * np.conj(b[idx2]) * c
    return float(np.real(k))
    
def finesse_bullseye_photodiode(maxtem, w=1.0, r=0.5887050112577, name='bullseye'):
    """Prints beat coefficients for HG modes on a bullseye photo detector
    in the format for insertion into the kat.ini file of Finesse.

    The default kat.ini numbers have been generated with:
    finesse_bullseye_photodiode(6)

    maxtem: highest mode order (int)
    w: beam radius on diode [m]
    r: radius of inner disk element [m]
    name: name of entry in kat.ini (string)

    The standard parameter of w=1 and r=0.5887050112577 have been chosen to achieve
    a small a HG00xHG00 coefficient of <1e-13.
    """
    assert (isinstance(maxtem, int) and maxtem >=0), "maxtem must be integer >=0" 

    print('PDTYPE {0}'.format(name))

    for n1 in np.arange(0,maxtem+1):
        for m1 in np.arange(0,maxtem+1):
            for n2 in np.arange(n1,maxtem+1):
                for m2 in np.arange(m1,maxtem+1):
                    c = HG_bullseye_coefficients_numerical(n1,m1, n2, m2, w, r)
                    if (np.abs(c)>1e-13 and not np.isnan(c)):
                        print("{:2d} {:2d} {:2d} {:2d} {:.15g}".format(n1,m1,n2,m2, c))
    print('END')
