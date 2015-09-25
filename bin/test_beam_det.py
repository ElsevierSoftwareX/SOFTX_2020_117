import pykat

kat = pykat.finesse.kat()

kat.parseCommands("""
l l1 1 0 n0
s s1 1 n0 n1
beam b1 0 n1

gauss g1 l1 n0 1 0

xaxis b1 x lin -3 3 100
x2axis b1 y lin -3 3 100
""")

out = kat.run()

import pylab
pylab.pcolormesh(out.x, out.y, out["b1"])

pylab.show()

