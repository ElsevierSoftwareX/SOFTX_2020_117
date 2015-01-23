import pykat

kat = pykat.finesse.kat()

kat.parseCommands("""
l l1 1 0 n0
s s1 100 n0 n1
m m1 1 0 0 n1 n2
s s2 200 n2 n3
bs bs1 1 0 45 45 n3 dump n4 dump
s s3 200 n4 n5
m m2 1 0 0 n5 n6
""")

kat.optivis().show()

