from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from pykat import finesse
from pykat.detectors import *
from pykat.components import *
from pykat.commands import *
from pykat.structs import *
import pylab as pl
import numpy as np
from pykat.parallel import parakat

# file uses parakat, start cluster with
# ipcluster-3.4 start -n=4
# or similar

# xbeta or ybeta
dof = "ybeta"

# switch between mirror and beamsplitter cavity
use_bs = True

ybeta = np.logspace(-8,-5, 100)

kat = finesse.kat()
kat.verbose = False
kat.maxtem = 4
kat.noxaxis = True

code_m = """
const fm 9M
l psl 1.0 0 npsl 
mod EOM $fm 0.001 1 pm 0 npsl nEOM1
s s1 0 nEOM1 nITM1

m1 ITM 0.02 0 0 nITM1 nITM2	 
attr ITM Rc -2500 
s s_cav 5000 nITM2 nETM1	
m1 ETM 0 0 0 nETM1 dump
attr ETM Rc 2700
cav c1 ITM nITM2 ETM nETM1

pd1 pdh1 $fm 0 nITM1
pd2 tf1 $fm 0 $fs nITM1
deriv_h 1e-13
"""

code_bs = """
const fm 9M
l psl 1.0 0 npsl 
mod EOM $fm 0.001 1 pm 0 npsl nEOM1
s s1 0 nEOM1 nITM1

# Cavity
# -------
bs1 ITM 0.02 0 0 0 nITM1 nITM1r nITM2 nITM2b  
attr ITM Rc -2500
s s_cav_a 5000 nITM2 nETM1
s s_cav_b 5000 nITM2b nETM1b	
bs1 ETM 0 0 0 10 nETM1 nETM1b dump dump
attr ETM Rc 2700 
cav c1 ITM nITM2 ITM nITM2b

# Output
# --------
pd1 pdh1 $fm 0 nITM1r
pd2 tf1 $fm 0 $fs nITM1r
deriv_h 1e-13
"""

if use_bs:
	kat.parseCommands(code_bs)
else:
	kat.parseCommands(code_m)

	
tf = np.zeros(len(ybeta))
pdh = np.zeros(len(ybeta))

pk = parakat()

for k, beta in enumerate(ybeta):
	setattr(kat.ETM, dof, beta)
	
	k1 = kat.deepcopy()
	k1.parseCommands('diff ETM %s' % dof)
	pk.run(k1, cmd_args=["-cr=on"])
	
	k2 = kat.deepcopy()
	k2.parseCommands('fsig sig1 ETM %s 1e-3 0 1' % dof)
	pk.run(k2, cmd_args=["-cr=on"])
	
outs = pk.getResults()
pk.close()

for k, beta in enumerate(ybeta):
	out1, out2 = outs[2*k:(2*k+2)]
	tf[k] = out2['tf1']
	pdh[k] = out1['pdh1']
	
# Plotting
# --------------------------------------------------

fig = pl.figure()

pl.loglog(ybeta, np.abs(pdh), label='Slope')
pl.loglog(ybeta, np.abs(tf),  label='TF(DC)')

pl.xlabel('Misalignment offset [rad] (divergence ~20urad)')
pl.ylabel('Error signal slope [W/rad]')
pl.legend(loc=0)

pl.title("%s - %s: signal frequency = %g Hz, maxtem = %i" % (("Beamsplitter" if use_bs else "Mirror"), dof, k2.signals.f.value, kat.maxtem))
pl.tight_layout()
pl.grid()
pl.show()
