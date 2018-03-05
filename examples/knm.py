import pykat

from pykat.optics.knm import plot_knm_matrix, knmHG, makeCouplingMatrix
from pykat.optics.maps import surfacemap

q1 = pykat.BeamParam(w0=5e-2, z=0)
q2 = pykat.BeamParam(w0=5e-2, z=0)

R  = max(q1.w, q2.w)

N  = 500
dx = R*6/N

m = surfacemap("empty",
               "phase both",
               size=(N, N),
               center=((N+1)/2, (N+1)/2),
               step_size=(dx, dx))

C = makeCouplingMatrix(2)
K = knmHG(C, q1, q2, surface_map=m)

plot_knm_matrix(C, abs(K))



C = [0,0,0,0] # 00 -> 00
K = knmHG(C, q1, q2, surface_map=m)

print(K)


C = [[0,0,1,0], [0,0,0,0]] # [00->10, 00->00]
K = knmHG(C, q1, q2, surface_map=m)

print(K)
