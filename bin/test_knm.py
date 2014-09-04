from pykat import *
from pykat.utilities.knm import knmHG, makeCouplingMatrix, plot_knm_matrix
from pykat.utilities.maps import aperturemap
import numpy as np
import time

q1 = beam_param(w0=5e-2, z=0)
q2 = beam_param(w0=5e-2, z=0)

aperture = 0.1
s = 1000
size = np.array([s, s])
stepsize = 0.3/(size-1)
smap = aperturemap("tilt", size, stepsize, aperture)

couplings = makeCouplingMatrix(1)

params = {"usepolar":True, "aperture":aperture, "epsabs": 1e-3, "epsrel": 1e-3}

t0 = time.time()
kbh = knmHG(couplings, q1, q2, method="adaptive", verbose=True, params=params)
print time.time() - t0

t0 = time.time()
kr = knmHG(couplings, q1, q2, surface_map=smap, method="riemann", verbose=True)
print time.time() - t0

smap.generateROMWeights(isModeMatched=True, verbose=True)
t0 = time.time()
krm = knmHG(couplings, q1, q2, surface_map=smap, method="romhom")
print time.time() - t0

print kbh
print kr
print krm

plot_knm_matrix(couplings, np.log10(np.abs(kbh - kr)))

print params