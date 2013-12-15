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

kat.add(photodiode('pd_ref','n3', num_demods=1, demods=[10,0]))

kat.add(xaxis("lin", [0, 1000], "pd_ref", "f1", 100))

out = kat.run(printout=0, printerr=0)

pl.figure()
pl.plot(out.x, out["pd_ref"])
pl.xlabel(out.xlabel)
pl.ylabel("Intensity [W]")
pl.legend(out.ylabels)
pl.show()
