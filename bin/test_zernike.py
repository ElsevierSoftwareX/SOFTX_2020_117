import numpy as np
from pykat.utilities.maps import zernikemap
import random
import pykat
import pylab
import os

r = 0.0125
size = np.array([200, 200])

pylab.figure()
legends = []
for n in xrange(0,3):
	for m in range(-n, n+1):
		legends.append("%i,%i" % (m,n))
		
		zmap = zernikemap('dsad', size, 2*r / size, 2*r, 1e-7)
		zmap.setZernike(m, n, 10)		
		zmap.type = 'phase reflection'
		zmap.write_map("zernike.map")

		kat = pykat.finesse.kat()

		kat.parseCommands("""
		l l1 1 0 n0
		m m1 0.99 0.01 0 n0 n1
		s s1 12 n1 n2
		m m2 0.999 0.001 0 n2 n3
		attr m1 Rc -20
		attr m2 Rc 20

		cav c1 m1 n1 m2 n2

		pd refl n0
		pd circ n1
		pd tran n3
		bp wz x w n2
		trace 2

		yaxis re:im
		map m2 %s/zernike.map

		maxtem 4
		xaxis m2 phi lin -10 190 2000
		""" %  os.path.realpath(os.path.curdir))

		out = kat.run(printerr=1)

		print out["wz"][0]
		
		pylab.semilogy(out.x, out["tran"])

pylab.legend(legends)
pylab.xlim(min(out.x), max(out.x))
pylab.show()
