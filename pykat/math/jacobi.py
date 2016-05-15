from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import scipy.special

def jacobi(n,a,b,x):
  """
  Jacobi Polynomial P_n^{a,b} (x)for real x.

                  n    / n+a \ / n+b \ / x-1 \^(n-s) / x+1 \^s
  P_n^{a,b}(x)= Sum    |     | |     | | --- |       | --- |
                  s=0   \ s  / \ n-s / \  2  /       \  2  /
  
  n,a,b (int)
  x (real)
  P (real)

  Implementation of Jacobi fucntion using binominal coefficients.
  This can handle values of alpha, beta < -1 which the special.eval_jacobi
  function does not.
  Andreas Freise 15.05.2016
  """

  P=0.0
  for s in np.arange(0,n+1):
    P=P+scipy.special.binom(n+a,s) * scipy.special.binom(n+b,n-s) * (x-1.0)**(n-s) * (x+1.0)**s
  P=P*0.5**n
  return P
