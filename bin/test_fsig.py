from pykat import finesse
from pykat.detectors import *
from pykat.components import *
from pykat.commands import *
from pykat.structs import *

import numpy as np
import pylab as pl

code = """
l l1 2 0 n1
m m1 0.99 0.01 0 n1 n2
s cav1 1200 n2 n3
m m2 0.99 0.01 -0.1 n3 n4

attr m2 m 1  # mech sus1

ad up_refl 0 n1
ad low_refl 0 n1

qd refl_A 0 0 n1
qd refl_Q 0 90 n1
qd tran_A 0 0 n4
qd tran_Q 0 90 n4

put up_refl f $x1
put low_refl f $mx1

yaxis log re:im

fsig noise m2 1 0
"""

kat = finesse.kat(kat_code=code)

kat.removeLine("fsig noise 9")

kat.signals.apply(kat.l1.P, 1, 0)
kat.signals.apply(kat.m1.phi, 1, 90)

kat.add(xaxis('log', [1, 1000], kat.signals.f, 100))

out = kat.run(printout=0, printerr=0)

# using real and imag part compute the complex value of the upper and lower sidebands
a_up = out.y[:,0] + out.y[:,1]*1j
a_lo = out.y[:,2] + out.y[:,3]*-1j

pl.figure(1)
pl.loglog(out.x, np.abs(a_up + a_lo), out.x, np.abs((a_up - a_lo) / (1j)))
pl.xlabel(out.xlabel)
pl.title("Reflection quadratures with no relative carrier phase")
pl.legend(["Amplitude","Phase"])
#pl.show()

kat.remove(kat.signals)