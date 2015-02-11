import pykat

kat = pykat.finesse.kat()

kat.parseCommands("""
l l1 1 0 n0
s s1 100 n0 n1
m m1 1 0 0 n1 n2
s s1 50 n2 n3

bs bs1 1 0 45 45 n3 n4 n5 n6

s s2 200 n5 n5a
m m2 1 0 0 n5a dump

s s3 200 n4 n4a
m m3 1 0 0 n4a dump

s s4 50 n6 n6a
m m4 1 0 0 n6a n6b
""")

kat.optivis().show()

