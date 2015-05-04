from pykat import finesse
from pykat.detectors import *
from pykat.components import *
from pykat.commands import *
from pykat.structs import *
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d.axes3d import Axes3D

import numpy as np
import pylab as pl

code = """
l l1 1 0 0 n1
s s1 10 1 n1 n2

ad ad1 0 n2
beam b1 0 n2
maxtem 0 

xaxis b1 x lin -10.0 10 100
x2axis b1 y lin -6 6 100

"""

kat = finesse.kat()

kat.parseCommands(code)
kat.s1.n1.q = pykat.beam_param(w0=1e-3, z=0)

out = kat.run(printout=0,printerr=0)

fig = pl.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
x, y = np.meshgrid(out.x, out.y)
p = ax.plot_wireframe(x, y, out["b1"])   
pl.xlabel(out.xlabel)
pl.ylabel(out.ylabel)
pl.show()

pl.figure()
extent = [np.min(out.x),np.max(out.x),np.min(out.y),np.max(out.y)]
imgplot = pl.imshow(out["b1"], extent=extent)
#imgplot.set_interpolation('bicubic')
imgplot.set_interpolation('nearest')
pl.colorbar()
pl.xlabel(out.xlabel)
pl.ylabel(out.ylabel)
pl.show()
