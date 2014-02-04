from pykat import finesse
from pykat.detectors import *
from pykat.components import *
from pykat.commands import *
from pykat.structs import *

import numpy as np
import pylab as pl

code = """
l l1 1 0 0 n1 ### test
s s1 10 1 n1 n2
m m1 0.5 0.5 0 n2 n3
s s2 10 1 n3 n4
m m2 0.5 0.5 0 n4 n5
s s3 10 1 n5 n6

yaxis abs:deg
"""

kat = finesse.kat()

kat.parseCommands(code)

kat.add(cavity('cav1', 'm1', 'n3', 'm2', 'n4'))

kat.add(photodiode('pd_ref','n2'))
kat.add(photodiode('pd_trs','n5'))
kat.add(photodiode('pd_cav','n4', num_demods=1, demods=[1]))


kat.add(xaxis("lin", [0, 360], kat.m2.phi, 100))

kat.m1.Rcx = -1000.0
kat.m1.Rcy = -1000.0
kat.m2.Rcx =  1000.0
kat.m2.Rcy =  1000.0

kat.maxtem = 0

out = kat.run(printout=0,printerr=0)

pl.figure()
pl.plot(out.x, out["pd_cav"])
pl.xlabel(out.xlabel)
pl.ylabel("Intensity [W]")
pl.legend(out.ylabels)
pl.show()
