from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import scipy.special
from math import factorial

def laguerre(p,l,x):
    """ function to evaluate associated Laguerre Polynomial L_p^l (x).
    Usage: L = laguerre (p,l,x)

                   p    1   / l+p \       
    L_p^l(x)=   Sum    ---  |     |  (-x)^j
                 j=0   j!   \ p-j /    
    p,l (int)
    x (real)
    L (real)
    Andreas Freise 15.05.2016
    """

    L=0.0
    for j in np.arange(0,p+1):
        L = L + scipy.special.binom(l+p,p-j) / factorial(j) * (-x)**j
    
    return L

