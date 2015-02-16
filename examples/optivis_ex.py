import pykat

kat = pykat.finesse.kat()

kat.parseCommands("""
l l1 1 0 n0
s s1 100 n0 n3

bs bs1 1 0 0 45 n3 n4 n5 n6

s s3 50 n4 n4a
m m1 1 0 0 n4a n4b
s s3a 200 n4b n7a
m m2 1 0 0 n7a n7b

s s4 50 n5 n5a
m m3 1 0 0 n5a n5b
s s4a 200 n5b n8a
m m4 1 0 0 n8a n8b
""")

kat.optivis().show()

