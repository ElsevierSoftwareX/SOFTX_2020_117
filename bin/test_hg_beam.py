import pykat
from pykat.utilities.optics.gaussian_beams import *

gx = beam_param(w0=2e-3, z=0)
gy = beam_param(w0=1e-3, z=0)

beam = HG_beam(gx,gy,0,0)
beam.n = 5
beam.m = 6

beam.plot()


