import pykat
import pylab

kat = pykat.finesse.kat()

kat.parseCommands("""
l l1 1 0 n0
m m1 0.9 0.1 0 n0 n1
s s1 1 n1 n2
m m2 0.9 0.1 0 n2 n3
attr m1 Rc -2
attr m2 Rc 2
pd circ n2
xaxis m1 phi lin 0 180 1000
yaxis log abs
maxtem 0
cav c1 m1 n1 m2 n2
""")

kat.timeCode = True

out, timings = kat.run()
