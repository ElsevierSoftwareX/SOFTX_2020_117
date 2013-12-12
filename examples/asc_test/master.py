from pykat import finesse
from pykat.commands import *
import pylab as pl

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





