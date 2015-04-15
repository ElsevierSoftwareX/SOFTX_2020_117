"""
Test file to ensure that references between deep copied objects are handled properly
"""

import pykat
from copy import deepcopy

kat0 = pykat.finesse.kat()

kat0.parseCommands("m m1 1 0 0 n0 n1")
kat0.parseCommands("pd o0 n1")
kat0.parseCommands("pd1 o1 1 0 n1")

kat1 = deepcopy(kat0)

assert(kat0 != kat1)
assert(kat0.__class__ != kat1.__class__)

assert(kat0.m1 != kat1.m1)
assert(kat0.m1.__class__ != kat1.m1.__class__)

assert(kat0.m1.n0 != kat1.m1.n0)
assert(kat0.m1.n0.node != kat1.m1.n0.node)
assert(kat0.nodes.n0 != kat1.nodes.n0)

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