import pykat
from pykat.tools.plotting.beamtrace import plot_beam_trace
import numpy as np
import pylab
import copy

kat = pykat.finesse.kat()

cmds = """
l l1 1 0 n0
s s0 100 n0 n1
m m1 0.5 0.5 0 n1 n2
s s1 100 n2 n3
m m2 0.5 0.5 0 n3 n4
s s2 20 n4 n5
bs bs1 0.5 0.5 0 0 n5 n6 n7 n8

s s3 20 n6 n9

s s4 20 n7 n10

s s5 20 n8 n11

bs bs2 0.5 0.5 0 0 n9 n12 n13 n14

s s6 3 n12 n15

lens lens1 200 n15 n16

s s7 500 n16 n17

#gouy g1 x s0 s1 s2

#bp bp1 x w0 n5
gauss g1 l1 n0 8e-3 -1000 4e-3 -1200
noxaxis
maxtem 0 
"""

kat.parseCommands(cmds)

plot_beam_trace(kat, 'n0', 'n17')