# easy test to make sure beam tracing is tracking the gouy phase accumulation
# and shape along a path correctly. Here I just expand and lens a beam back again
# Gouy phase should be about 360 deg in total.
from __future__ import print_function

import numpy as np
import pykat

kat = pykat.finesse.kat()
kat.parse("""
l l1 1 0 n0
s s1 1000 n0 n1
bs m1 1 0 0 0 n1 n2 dump dump
s s2 1000 n2 n3
bs m2 1 0 0 0 n3 n4 dump dump
s s3 2000 n4 n5
""")

bp = pykat.BeamParam(w0=1e-5, z=-1000)

T = kat.beamTrace(bp, "n0", "n5")

kat.m2.Rc = T.data['m2']['qin'].Rc

T = kat.beamTrace(bp, "n0", "n5")

assert(np.allclose(359.99993233082711, T.data['s3']['gouy']))
assert(np.allclose(999.9999999999994+0.00029526246744265j, T.data['s3']['qout'].q))
