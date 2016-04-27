"""
--------------------------------------------------------------
Example file for using PyKat to remove objects from a 'kat'
object for Finesse simulations
Finesse: http://www.gwoptics.org/finesse
PyKat:   http://www.gwoptics.org/pykat

Example showing how to remove and add components with Pykat.
Commands that get parsed into a pykat object can be interacted 
with in a object orientated manner. So you can call kat.component.remove()
to remove it for example.

This example shows how a kat object can be manipulated when pykat objects
are available.

If the a Finesse command hasn't been implemented as a pykat object yet
the command will be added as an "extra line", these can still be removed
and added to using the kat.removeLine and kat.addLine commands.
You can set
  kat.verbose = True
to print a warning message when "extra lines" are generated
Daniel Brown 17/12/14
--------------------------------------------------------------
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pykat
from pykat.components import *

kat = pykat.finesse.kat()

kat.parseCommands("""
m m1 0.9 0.1 0 n1 n2
s s1 1 n2 n3
m m2 0.9 0.1 0 n3 n4
""")

print ("Before...")
print ("".join(kat.generateKatScript()))

kat.s1.remove()

print ("After remove...")
print ("".join(kat.generateKatScript()))

# Adding in with commands
kat.parseCommands("""
s s2 1 n2 n2a
m m3 0.9 0.1 0 n2a n3a
s s3 1 n3a n3
""")

print ("After add with commands...")
print ("".join(kat.generateKatScript()))

kat.s2.remove()
kat.s3.remove()
kat.m3.remove()

# Adding in with objects
kat.add(space("s2", "n2", "n2a", L=1))
kat.add(space("s3", "n3", "n3a", L=1))
kat.add(mirror("m3", "n2a", "n3a", R=0.9, T=0.1, phi=0))

print ("After add with objects...")
print ("".join(kat.generateKatScript()))
