# This file tests the bayer helms mismatch calculation against a numerical
# Newton-Cotes integration of the coupling coefficients

import pykat
from pykat.optics.knm import *
from pykat.optics.maps import *
import numpy as np

N = 501

dx = 1/float(N)

m = curvedmap("test", (N,N), (dx,dx), 1e15)

C = makeCouplingMatrix(0)

q1 = pykat.BeamParam(w0=5e-2, z=0)
q2 = pykat.BeamParam(w0=10e-2, z=10)

Kmap = knmHG(C, q1, q2, surface_map=m, method="riemann", cache=False)
Kbh  = knmHG(C, q1, q2, method="bayerhelms", cache=False)

print(Kmap)
print(Kbh)

print(np.max(abs(Kmap - Kbh)))

assert(np.max(abs(Kmap - Kbh)) < 1e-12)
