from pykat import finesse
from pykat.detectors import *
from pykat.components import *
from pykat.commands import *
from pykat.structs import *

import numpy as np
import pylab as pl

code = """
l l1 1 0 n1
m m1 0.99 0.01 0 n1 n2
s cav1 1200 n2 n3
m m2 0.99 0.01 -0.1 n3 n4

#attr m1 m 1# mech sus1
attr m2 m 1# mech sus1

fsig sig l1 amp 1 0 4

#ad car_refl 0 n1
ad up_refl 0 n1
ad low_refl 0 n1

#qd refl_A 0 0 n1
#qd refl_Q 0 90 n1
#qd tran_A 0 0 n4
#qd tran_Q 0 90 n4
put up_refl f $x1
put low_refl f $mx1

xaxis sig f log 1 10000 1000
yaxis log re:im
"""

kat = finesse.kat(kat_code=code)
run = kat.run(printout=0,printerr=0)

# using real and imag part compute the complex value of the upper and lower sidebands
a_up = run.y[:,0] + run.y[:,1]*1j
a_lo = run.y[:,2] + run.y[:,3]*-1j

pl.figure()
pl.loglog(run.x, np.abs(a_up + a_lo), run.x, np.abs((a_up - a_lo) / (1j)))
pl.xlabel(run.xlabel)
pl.show()

pl.figure()
pl.loglog(run.x, np.abs(a_up), run.x, np.abs(a_lo))
pl.xlabel(run.xlabel)
pl.show()