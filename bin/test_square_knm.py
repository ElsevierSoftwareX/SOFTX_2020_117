import pykat
import numpy as np
from pykat.optics.knm import square_aperture_HG_knm, riemann_HG_knm
import math

R = 0.15
q = pykat.beam_param(w0=0.05, z=0)
N = 1200

mode_i = (1,1)
mode_o = (1,1)

k = square_aperture_HG_knm(mode_i, mode_o, q, R)

x = y = np.linspace(-R, R, N)
Axy = np.ones((N,N))
k_ = riemann_HG_knm(x, y, mode_i, mode_o, q, q, Axy=Axy, newtonCotesOrder=1)

print "%15.15f + %15.15fi" % (k.real, k.imag)
print "%15.15f + %15.15fi" % (k_.real, k_.imag), abs(k-k_)/1e-6, "ppm"