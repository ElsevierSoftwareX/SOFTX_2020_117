from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import scipy.special

def jacobi(n,a,b,x):
  """ Implementation of Jacobi fucntion using binominal coefficients.
  This can handle values of alpha, beta < -1 which the special.eval_jacobi
  function cannot."""
  P=0.0
  for s in np.arange(0,n+1):
    P=P+scipy.special.binom(n+a,s) * scipy.special.binom(n+b,n-s) * (x-1.0)**(n-s) * (x+1.0)**s
  P=P*0.5**n
  return P
