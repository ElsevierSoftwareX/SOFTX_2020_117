from pykat.utilities.maps import read_map
from pykat.utilities.knm import *
from pykat.utilities.optics.gaussian_beams import beam_param

couplings = np.array([[[0,0,0,0], [0,0,1,0], [0,0,0,1]],
                      [[1,0,0,0], [1,0,1,0], [1,0,0,1]],
                      [[0,1,0,0], [0,1,1,0], [0,1,0,1]]])

q1 = beam_param(w0=2e-2, z=0)
q2 = beam_param(w0=2e-2, z=0)

m = read_map('etm08_virtual.txt')

m.data = np.zeros((200, 200))
m.x = np.linspace(-0.151, 0.149, m.size[0])
m.y = np.linspace(-0.151, 0.149, m.size[1])
m.center = np.array(m.size)/2
m.step_size = (m.x[1]-m.x[0], m.y[1]-m.y[0])

print "generating weights..."
w = m.generateROMWeights()

print "computing Riemann knm..."
K1 = knmHG(couplings, q1, q2, surface_map=m)
print np.abs(K1)

print "computing ROMHOM knm..."
K2 = knmHG(couplings, q1, q2, surface_map=m, method="romhom")
print np.abs(K2)

#w.writeToFile("testWeights.rom")




