import pykat

kat = pykat.finesse.kat()

kat.parseCommands("""
l l1 1 0 80 n1
l l2 1 0 90 n4

bs bs1 0.5 0.5 0 0 n1 n2 n3 n4

fsig noise l1 amp 1 0 1

hd hd1 0 n2 n3

xaxis hd1 phase lin 0 360 1000
""")