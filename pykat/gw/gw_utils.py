"""
------------------------------------------------------
Collection of various utility functions related to
designing interferometers for graviational wave
detection

Andreas 23.06.2016
http://www.gwoptics.org/pykat/
------------------------------------------------------
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np

def SQL_x(f_vector, m, N=2):
    """
    Function to compute the Standard Quantum Limit (SQL)
    as amplitude spectral density for mirror displacement.
    
    f_vector = array of frequencies (x-axis) [Hz]
    m = mass of test masses (mirrors) [kg]
    N = 1: Michelson without arm cavities
    N = 2: Michelson with arm cavities (default)

    See Kimble et al. `Conversion of conventional gravitational-wave
    interferometers into quantum nondemolition interferometers by
    modifying their input and/or output optics' Physical Review D,
    2002, 65, 022002-+
    """

    hbar=6.62606957E-34/(2.0 *np.pi)
    SQLx= np.sqrt( N * 4 * hbar / ( m * (2 * np.pi * f_vector)**2))
    return SQLx

def SQL_h(f_vector, m, Larm, N=2):
    """
    Function to compute the Standard Quantum Limit (SQL)
    as amplitude spectral density for GW strain.
    
    f_vector = array of frequencies (x-axis) [Hz]
    m = mass of test masses (mirrors) [kg]
    Larm = length of interferometer arms [m]
    N = 1: Michelson without arm cavities
    N = 2: Michelson with arm cavities (default)

    See Kimble et al. `Conversion of conventional gravitational-wave
    interferometers into quantum nondemolition interferometers by
    modifying their input and/or output optics' Physical Review D,
    2002, 65, 022002-+
    """
    
    SQLx=SQL_x(f_vector, m, N=N)
    return SQLx/Larm
