from pykat import finesse
from pykat.detectors import *
from pykat.components import *
from pykat.commands import *
from pykat.structs import *
from pykat.plotting import *
import numpy as np
import pylab as pl

kat = finesse.kat()

laser(kat,'l1','n1',1)
space(kat,'s1','n1','n2',1)

mirror(kat,'m1','n2','n3',R=0.8,T=0.2)
space(kat,'s2','n3','n4',L=1)
mirror(kat,'m2','n4','n5',R=0.7,T=0.3)
cavity(kat, 'cav1','m1','n3','m2','n4')
space(kat,'s3','n5','n6',L=1)

photodiode(kat,'pd_cav','n4')
photodiode(kat,'pd_ref','n2')
photodiode(kat,'pd_trs','n5')

kat.m1.Rcx = -1000.0
kat.m1.Rcy = -1000.0
kat.m2.Rcx =  1000.0
kat.m2.Rcy =  1000.0

xaxis(kat, Scale.linear, [0,360], kat.m2, kat.m2.phi, 1000)

kat.maxtem = 0

run = kat.run(printout=0,printerr=0)

pl.figure()
pl.plot(run.x,run.y)
pl.xlabel(run.xlabel)
pl.ylabel("Intensity [W]")
pl.legend(run.ylabels)
#pl.show()

kat.m1.R = 0.5
kat.m1.T = 0.5
kat.pd_cav.enabled = False

run = kat.run(printout=0,printerr=0)

pl.figure()
pl.plot(run.x,run.y)
pl.xlabel(run.xlabel)
pl.ylabel("Intensity [W]")
pl.legend(run.ylabels)
#pl.show()

kat.openGUI()


