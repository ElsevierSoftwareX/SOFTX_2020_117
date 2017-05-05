import pykat
from pykat.optics.knm import *
from pykat.optics.maps import *


m = read_map("/Users/ddb/git/aligo_finesse/PI/20153006/maps/mode_37.map")
C = makeCouplingMatrix(10)

q1 = pykat.BeamParam(w0=5e-2, z=0)

K = knmHG(C, q1, q1, surface_map=m)
plot_knm_matrix(C, abs(K))