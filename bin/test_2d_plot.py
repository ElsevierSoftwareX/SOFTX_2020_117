from pykat import finesse
from pykat.utilities.optics.gaussian_beams import gauss_param
from pykat.detectors import *
from pykat.components import *
from pykat.commands import *
from pykat.structs import *
from mpl_toolkits.mplot3d.axes3d import Axes3D

import numpy as np
import pylab as pl

code = """
l l1 1 0 0 n1
s s1 10 1 n1 n2

beam b1 0 n2
maxtem 0 

xaxis b1 x lin -7.0 7.0 50
x2axis b1 y lin -7.0 7.0 50

"""

kat = finesse.kat()

kat.parseCommands(code)
kat.s1.n1.q = gauss_param(w0=1e-3, z=0)

out = kat.run(printout=0,printerr=0)

fig = pl.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
x, y = np.meshgrid(out.x, out.y)
p = ax.plot_surface(x, y, out.z[0,:])   
pl.show()
