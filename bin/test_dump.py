import pykat

kat = pykat.finesse.kat()

kat.parseCommands("""
bs bs1 0.5 0.5 0 0 n1 dump n3 dump
""")

print "BEFORE"
print "".join(kat.generateKatScript())

kat.nodes.replaceNode(kat.bs1, kat.bs1.nodes[3], kat.nodes.createNode("test4"))
kat.nodes.replaceNode(kat.bs1, kat.bs1.nodes[1], kat.nodes.createNode("test2"))

print "AFTER"
print "".join(kat.generateKatScript())
