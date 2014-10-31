import pykat

kat = pykat.finesse.kat()

kat.parseCommands("""
l l1 1 0 0 n1
sq l2 0 10 0 n4

bs bs1 0.5 0.5 0 0 n1 n2 n3 n4

fsig noise l1 amp 1 0 1

qhd qhd180 180 n2 n3
qhd qhd0 0 n2 n3

xaxis l1 phase lin 0 360 360
""")

out = kat.run()

#out.plot()

print kat.qhd180
