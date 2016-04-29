import pykat
from copy import deepcopy

kat = pykat.finesse.kat()

code_det = """
m m1 1 0 0 n0 n1
pd1 pdr 9M 90 n1
"""

kat.parseCommands(code_det)

kat.pdr.f1 = "0.1k"
assert(kat.pdr.f1 == 100)
assert(type(kat.pdr.f1) is pykat.param.Param)

kat.pdr.phase1 = "10u" 
assert(kat.pdr.phase1 == 1e-5)
assert(type(kat.pdr.phase1) is pykat.param.Param)

kat.m1.R = "10000u"
assert(kat.m1.R == 0.01)

#################################
kat = deepcopy(kat)

kat.m1.R = 0.9
assert(kat.m1.R == 0.9)

kat.pdr.phase1 = 20
assert(kat.pdr.phase1 == 20)

print("PASSED")