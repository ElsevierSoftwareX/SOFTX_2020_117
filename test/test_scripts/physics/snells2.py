# Testing Snells law
# This ones checks that the angles between mirrors and beamsplitters should cancel out
# and be computed correctly from HG10 amplitudes

import pykat

kat = pykat.finesse.kat()

kat.parse("""
l l1 1 0 nb
bs bs1 1 0 0 0 na nb dump dump

attr bs1 xbeta 0.5n
attr m1  xbeta -1n

s s0 0 1 na n0
m m1 0 1 0 n0 n1
s s1 0 2 n1 n2

gauss g1 s0 n0 0.1 0
maxtem 1

noxaxis

ad a10 1 0 0 n0*
ad b10 1 0 0 n2
""")

out, trace = kat.run(getTraceData=True)

qa = trace[0]['n0'][0]
qb = trace[0]['n2'][0]

diva = qa.w0/qa.zr
divb = qb.w0/qb.zr

tilt = lambda a, q: a * (q.w0/q.zr) # tilt is amplitude * divergence

tilta = tilt(out['a10'], qa)
tiltb = tilt(out['b10'], qb)

assert(abs((tilta + (kat.m1.xbeta.value-tilta) * kat.s0.n.value/kat.s1.n.value) - tiltb) < 1e-15)
