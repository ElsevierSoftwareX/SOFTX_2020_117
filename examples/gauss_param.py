from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from pykat.optics.ABCD import apply, mirror_trans
from pykat.optics.gaussian_beams import beam_param

nr1 = 1
nr2 = 1.44963098985906
q1 = beam_param(q=5.96343 + 3.04713j)
abcd = mirror_trans(nr1, nr2, float("inf"))
# into material
q2 = apply(abcd, q1.q, nr1, nr2)
# and out again
q3 = apply(abcd, q2.q, nr2, nr1)

print ("q1 =", q1, " w0 =", q1.w0, " w =", q1.w, " z =", q1.z)
print ("q2 =", q2, " w0 =", q2.w0, " w =", q2.w, " z =", q2.z)
print ("q3 =", q3, " w0 =", q3.w0, " w =", q3.w, " z =", q3.z)
