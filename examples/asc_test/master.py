from pykat import finesse
from pykat.commands import *
import pylab as pl
import copy

print """
--------------------------------------------------------------
Example file for using PyKat to automate Finesse simulations
Finesse: http://www.gwoptics.org/finesse
PyKat:   https://pypi.python.org/pypi/PyKat/

The file runs through the various pykat files which are used
to generate the Finesse results reported in the document:
`Comparing Finesse simulations, analytical solutions and OSCAR 
simulations of Fabry-Perot alignment signals', LIGO-T1300345

Andreas Freise 06.12.2013
--------------------------------------------------------------
"""


kat = finesse.kat(tempdir=".",tempname="test")
#kat = finesse.kat()
kat.verbose = False
kat.loadKatFile('asc_base.kat')
kat.maxtem=3
Lambda=1064.0e-9


print "--------------------------------------------------------"
print " 1. tunes ETM position to find resonance"
import asc_resonance 
kat.ETM.phi=asc_resonance.run(kat)

print "--------------------------------------------------------"
print " 2. print sideband and carrier powers/amplitudes"
import asc_powers 
asc_powers.run(kat)

print "--------------------------------------------------------"
print " 3. determine the optimal phase for the PDH signal"
import asc_pd_phase
(p_phase, q_phase) = asc_pd_phase.run(kat)

# setting demodulation phase
code_det = """
pd1 PDrefl_p 9M 0 nWFS1
scale 2 PDrefl_p
pd1 PDrefl_q 9M 90 nWFS1
scale 2 PDrefl_q
"""
kat.parseKatCode(code_det)
kat.PDrefl_p.phi[0]=p_phase
kat.PDrefl_q.phi[0]=q_phase

print "--------------------------------------------------------"
print " 4. adding a 0.1nm offset to ETM and compute PDH signal"
phi0=kat.ETM.phi
kat.ETM.phi=phi0 + 0.1/1064.0*360
print " new ETM phi tuning = %g " % kat.ETM.phi

import asc_pd_signal
(pd_p, pd_q) = asc_pd_signal.run(kat)
print " PDH inphase     = %e " % pd_p
print " PDH quadrtature = %e " % pd_q







