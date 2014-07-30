from pykat.utilities.maps import read_map
from pykat.utilities.knm import *
from pykat.utilities.optics.gaussian_beams import beam_param
import time

couplings = np.array([[[0,0,0,0], [0,0,1,0], [0,0,0,1]],
                      [[1,0,0,0], [1,0,1,0], [1,0,0,1]],
                      [[0,1,0,0], [0,1,1,0], [0,1,0,1]]])

q1 = beam_param(w0=4e-2, z=0)
q2 = beam_param(w0=4e-2, z=0)

m = read_map('etm08_virtual.txt')

m.data = np.zeros((100, 100))
m.x = np.linspace(-0.15, 0.15, m.size[0])
m.y = np.linspace(-0.15, 0.15, m.size[1])
m.center = np.array(m.size)/2
m.step_size = (m.x[1]-m.x[0], m.y[1]-m.y[0])

print "generating weights..."
t0 = time.time()
w, EI = m.generateROMWeights(isModeMatched=True)
print "Completed in ", time.time()-t0

print "computing Riemann knm..."
t0 = time.time()
K1 = knmHG(couplings, q1, q2, surface_map=m)
tr = time.time()-t0
print "Completed in ", tr
print np.abs(K1)

print "computing ROMHOM knm..."
t0 = time.time()
K2 = knmHG(couplings, q1, q2, surface_map=m, method="romhom")
tt = time.time()-t0
print "Completed in ", tt
print np.abs(K2)

print "Speed up", tr/tt

#w.writeToFile("testWeights.rom")




