import pykat
import numpy as np
from pykat.optics.knm import square_aperture_HG_knm, riemann_HG_knm
import math

R = 0.15
q = pykat.beam_param(w0=6e-3, z=2000)
N = 1201

mode_i = (0,0)
mode_o = (0,0)

k = square_aperture_HG_knm(mode_i, mode_o, q, R)

x = y = np.linspace(-R, R, N)
Axy = np.ones((N,N))
k_ = riemann_HG_knm(x, y, mode_i, mode_o, q, q, Axy=Axy, newtonCotesOrder=2)

print "%15.15f + %15.15fi" % (k.real, k.imag)
print "%15.15f + %15.15fi" % (k_.real, k_.imag), abs(k-k_)/1e-6, "ppm"