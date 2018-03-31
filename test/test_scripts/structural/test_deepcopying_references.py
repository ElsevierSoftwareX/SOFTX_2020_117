"""
Test file to ensure that references between deep copied objects are handled properly
"""

import pykat
from copy import deepcopy

kat0 = pykat.finesse.kat()

kat0.parse("m m1 1 0 0 n0 n1")
kat0.parse("pd1 o1 1 0 n1")

kat1 = deepcopy(kat0)

assert(kat0 != kat1)
assert(kat0.__class__ != kat1.__class__)

assert(kat0.m1 != kat1.m1)
assert(kat0.m1.__class__ != kat1.m1.__class__)

assert(kat0.m1.n0 != kat1.m1.n0)
assert(kat0.m1.n0.node != kat1.m1.n0.node)
assert(kat0.nodes.n0 != kat1.nodes.n0)

assert(kat0.o1 != kat1.o1)
assert(kat0.o1.__class__ != kat1.o1.__class__)

assert(kat0.m1.phi.owner == kat0.m1)
assert(kat1.m1.phi.owner == kat1.m1)

# use is to compare if two params are the same object as equals is override to compare the value
assert(kat0.o1.f1 is not kat1.o1.f1)
assert(kat0.o1.f1 == kat1.o1.f1)

kat0.o1.f1 *= 2
kat0.o1.phase1 *= 2

assert(isinstance(kat0.o1.f1, pykat.param.Param))
assert(isinstance(kat0.o1.phase1, pykat.param.Param))
assert(kat0.o1.f1 != kat1.o1.f1)

kat1.o1.num_demods = 2

assert(hasattr(kat1.o1, "f2"))
assert(not hasattr(kat0.o1, "f2"))

kat1.o1.num_demods = 1

assert(hasattr(kat1.o1, "f1"))

new = kat1.nodes.createNode("n4")

kat1.nodes.replaceNode(kat1.m1, kat1.m1.n0, new)

assert(not hasattr(kat1.nodes, "n0"))
assert(hasattr(kat1.nodes, "n4"))

assert(hasattr(kat0.nodes, "n0"))
assert(not hasattr(kat0.nodes, "n4"))

assert(hasattr(kat1.m1, "n4"))
assert(not hasattr(kat1.m1, "n0"))

assert(hasattr(kat0.m1, "n0"))
assert(not hasattr(kat0.m1, "n4"))

assert(kat0.nodes.n0 == kat0.m1.n0.node)
assert(kat0.nodes.n1 == kat0.m1.n1.node)

assert(kat1.nodes.n4 == kat1.m1.n4.node)
assert(kat1.nodes.n1 == kat1.m1.n1.node)

print("PASSED")