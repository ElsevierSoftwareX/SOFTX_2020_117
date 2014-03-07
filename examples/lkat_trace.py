import pykat
								 
cmd = """
l l1 1 0 n1
s s1 1 n1 n2
m m1 0.99 0.01 0 n2 n3
s s2 999 n3 n4
m m2 0.99 0.01 0 n4 n5
pd circ n3

noxaxis
maxtem 0

attr m1 Rc 0
attr m2 Rc 1000 
cav c1 m1 n3 m2 n4
"""

kat = pykat.finesse.kat()

kat.parseCommands(cmd)

info = kat.lkat_trace()

print "n1 qx =",info["n1"].qx

print "Cavity info ", info["c1"]

