#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pykat import finesse
from pykat.detectors import *
from pykat.components import *
from pykat.commands import *
from pykat.structs import *

import numpy as np
import pylab as pl

code = """
l l1 2 0 n1
m m1 0.99 0.01 0 n1 n2
s cav1 $test n2 n3
m m2 0.99 0.01 -0.1 n3 n4

attr m2 m 1  # mech sus1

const test 1200
ad up_refl 0 n1
ad low_refl 0 n1

qd refl_A 0 0 n1
qd refl_Q 0 90 n1
qd tran_A 0 0 n4
qd tran_Q 0 90 n4

put up_refl f $x1
put low_refl f $mx1

yaxis log re:im

fsig noise 1
"""

kat = finesse.kat(kat_code=code)

kat.signals.apply(kat.l1.P, 1, 0)
kat.signals.apply(kat.m1.phi, 1, 90)

kat.add(xaxis('log', [1, 1000], kat.signals.f, 100))

out = kat.run(printout=0, printerr=0)
