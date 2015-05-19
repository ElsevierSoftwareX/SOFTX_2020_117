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

ad up_refl $fs n1
ad low_refl $fs n1

qd refl_A 0 0 n1
qd refl_Q 0 90 n1
qd tran_A 0 0 n4
qd tran_Q 0 90 n4

qnoised qnd 1 $fs max n4
qshot   qsd 1 $fs max n4

pd1 p1 $fs max n4

fsig noise l1 amp 1 0
"""

kat = finesse.kat(kat_code=code)

# kat.signals.apply(kat.m2.z, 1.0, 0.0)
# kat.yaxis = "re:im"
# kat.add(xaxis('log', [1, 1000], kat.signals.f, 1000))
#
# out = kat.run(printout=0, printerr=0)
#
# a_up = out[kat.up_refl]
# a_lo = out[kat.low_refl]
#
# pl.figure(1)
# ax = pl.subplot(111)
# pl.grid()
# ax.loglog(out.x, abs(a_up))
# ax = ax.twinx()
# ax.plot(out.x, np.rad2deg(np.angle(a_up)), 'r')
# pl.xlabel(out.xlabel)

#pl.show()
print(kat.generateKatScript())
kat.remove(kat.signals)
