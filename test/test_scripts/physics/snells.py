# Testing Snells law

import pykat

kat = pykat.finesse.kat()

kat.parse("""
l l1 1 0 ni
tem l1 0 0 1 0
tem l1 1 0 1u 0

s s0 1 1 ni n0
m m1 0 1 0 n0 n1
s s1 1 2 n1 n2

gauss g1 s0 n0 0.1 0
maxtem 1

noxaxis

ad a10 1 0 0 n0*
ad b10 1 0 0 n1
""")

out, trace = kat.run(getTraceData=True)

qa = trace[0]['n0'][0]
qb = trace[0]['n1'][0]
print(qa.q, qb.q, qa.nr, qb.nr)

diva = qa.w0/qa.zr
divb = qb.w0/qb.zr
print(diva, divb)

print(out['a10'], out['b10'])

tilt = lambda a, q: a * (q.w0/q.zr) # tilt is amplitude * divergence

tilta = tilt(out['a10'], qa)
tiltb = tilt(out['b10'], qb)

print(tilta, tiltb, tilta/tiltb)

# Snells
# tiltb = n1/n2 tilta
print(tilta * kat.s0.n.value/kat.s1.n.value)

assert(abs(tilta * kat.s0.n.value/kat.s1.n.value- tiltb) < 1e-15)