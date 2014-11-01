import pykat

kat = pykat.finesse.kat()

kat.parseCommands("""
bs bs1 0.5 0.5 0 0 n1 dump n3 dump
""")

print "BEFORE"
print "".join(kat.generateKatScript())

kat.nodes.replaceNode(kat.bs1, kat.bs1.nodes[3], kat.nodes.createNode("test4"))
kat.nodes.replaceNode(kat.bs1, kat.bs1.nodes[1], kat.nodes.createNode("test2"))

kat.nodes.replaceNode(kat.bs1, "n1", kat.nodes.createNode("test1"))
kat.nodes.replaceNode(kat.bs1, "n3", kat.nodes.createNode("dump"))
kat.nodes.replaceNode(kat.bs1, "test1", kat.nodes.createNode("dump"))

print "AFTER"
print "".join(kat.generateKatScript())



