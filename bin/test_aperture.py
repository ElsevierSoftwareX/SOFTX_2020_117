from pykat import finesse
from pykat.commands import xaxis
import pylab as pl
import numpy as np

code = """
l l1 1 0 0 n1
s s1 10 1 n1 n2
m m1 1 0 0 n2 n3

gauss g1 m1 n2 1e-3 0

pd refl n2

xaxis m1 r_ap lin 0.1e-3 2e-3 10
"""

kat = finesse.kat(kat_code = code)

maxtem = np.arange(0, 1, 2)

for tem in maxtem:
    print "Calculating maxtem ",tem,"..."
    kat.maxtem = tem
    r = kat.run()
    pl.plot(r.x/1e-3, r.y, label="maxtem={0}".format(tem))

    
pl.ylabel("Reflected Power [W]")
pl.xlabel("Mirror aperture [mm]")
pl.legend()
pl.show()
    