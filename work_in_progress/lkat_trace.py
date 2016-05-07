from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import pykat
import pylab
import numpy
import ctypes
import pylibkat								 

cmd = """
l l1 1 0 n1
s s1 1 n1 n2
m m1 0.99 0.01 0 n2 n3
s s2 1 n3 n4
m m2 0.99 0.01 0 n4 n5
pd circ n3

noxaxis
maxtem 0

attr m1 Rc 0
attr m2 Rc 1000 
cav c1 m1 n3 m2 n4
retrace force

"""

kat = pykat.finesse.kat()

kat.parseCommands(cmd)

info = kat.lkat_trace()

print ("n1 qx =", info["n1"].qx)

print ("Cavity info ", info["c1"])

import numpy

Lengths = numpy.linspace(1, 999, 10000)
g_factors = []
        
lkat = ctypes.PyDLL("libkat.dylib")

try:
	lkat._pykat_preInit() # must always be called, sets up
	        # exception handling and such no simulation
	        # specifc code here

	# reads in the kat.ini and setups up other parts
	lkat._pykat_init()
	
	cmd = "".join(kat.generateKatScript())
	
	lkat._pykat_setup(cmd)
	
	inter = pylibkat.interferometer.in_dll(lkat, "inter")
	s2 = inter.space_list[1]
	cav = inter.cavity_list[0]
	
	inter.rebuild = 2
	inter.retrace = 1
	s2.rebuild = 2
	
	for L in Lengths:
		s2.L = L
		lkat._pykat_step()
		
		g_factors.append(cav.stability_x)

except Exception as ex: 
    print ("Exception caught in python: ", ex.message)
finally:
    # This should always be called no matter what
    lkat._pykat_finish(0)
	
pylab.plot(Lengths, g_factors)
pylab.show()


