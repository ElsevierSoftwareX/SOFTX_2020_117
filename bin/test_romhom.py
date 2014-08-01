from pykat.utilities.maps import read_map, tiltmap, surfacemap
from pykat.utilities.knm import *
from pykat.utilities.optics.gaussian_beams import beam_param
import time
import numpy as np

couplings = makeCouplingMatrix(5)

q1 = beam_param(w0=3e-2, z=0)
q2 = beam_param(w0=3e-2, z=0)

size = np.array([1200, 1200])
stepsize = 0.3/(size-1)

# This map has points evenly spaced about the axes
m = tiltmap("tilt", size, stepsize, (1e-6, 1e-6))
# Shifting the central point changes the results slightly
# but this is about a 1e-7 change, probably due to extra
# clipping
m.center += 0

print "Generating weights..."
t0 = time.time()
w, EI = m.generateROMWeights(isModeMatched=True)
w.writeToFile("testWeights.rom")
print "Completed in ", time.time()-t0

print "computing Riemann knm..."
t0 = time.time()
K1 = knmHG(couplings, q1, q2, surface_map=m)
tr = time.time()-t0
print "Completed in ", tr
#print np.abs(K1)

print "computing ROMHOM knm..."
t0 = time.time()
K2 = knmHG(couplings, q1, q2, surface_map=m, method="romhom")
tt = time.time()-t0
print "Completed in ", tt
#print np.abs(K2)
print "Speed up", tr/tt






