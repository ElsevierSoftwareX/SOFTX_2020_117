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
    
