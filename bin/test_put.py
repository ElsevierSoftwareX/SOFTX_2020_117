from pykat import finesse
from pykat.detectors import *
from pykat.components import *
from pykat.commands import *
from pykat.structs import *

import numpy as np
import pylab as pl

code = """
l l1 1 0 0 n1
s s1 10 1 n1 n2
mod eom 10 0.1 1 am 0 n2 n3
"""

kat = finesse.kat()

kat.parseCommands(code)

kat.add(pd('pdp',1,'n3'))
kat.add(pd('pdm',1,'n3'))

kat.add(xaxis("lin", [0, 1000], kat.eom, "f", 100))

kat.pdp.f1.put(kat.xaxis.x)
kat.pdm.f1.put(kat.xaxis.mx)

out = kat.run(printout=0, printerr=0)

pl.figure()
pl.plot(out.x, out["pdp"], out.x, out["pdm"])
pl.xlabel(out.xlabel)
pl.ylabel("Intensity [W]")
pl.legend(out.ylabels)
pl.show()
