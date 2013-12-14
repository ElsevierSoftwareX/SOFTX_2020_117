# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 15:29:19 2013

Takes a set of Finesse commands as input, parses it and outputs it again.
Used to check whether Finesse components, detectors, etc. are properly
reproducing the parameters they are given.

@author: Sean Leavey
"""

from pykat import finesse

code = """
l l1 1 0 0 n1
s s1 10 1 n1 n2
m m1 1 0 0 n2 n3
gr4 grating 1500 n4 n5 n6 n7
isol isolator 60 n8 n9
lens lens 10 n10 n11

pd refl n2

xaxis m1 r_ap lin 0.1e-3 2e-3 10
"""

kat = finesse.kat()
kat.parseCommands(code)

scriptList = kat.generateKatScript()

print ''.join(scriptList)