from pykat import finesse
from pykat.commands import xaxis
import pylab as pl

code = """
l l1 1 0 0 n1
s s1 10 1 n1 n2
m m1 0.5 0.5 0 n2 n3
s s2 10 1 n3 n4
m m2 0.5 0.5 0 n4 dump

ad ad1 0 n2
"""

kat = finesse.kat(kat_code = code)

kat.add(xaxis("lin", [0, 360], kat.m2, kat.m2.phi, 1000))

r = kat.run(printerr=1)

pl.plot(r.x, r.y)
pl.show()
    