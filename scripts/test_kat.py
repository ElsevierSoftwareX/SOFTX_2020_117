import sys
sys.path.append('../')

from pykat import finesse
from pykat.detectors import *
from pykat.components import *
from pykat.commands import *
from pykat.structs import *
from pykat.plotting import *

kat = finesse.kat()

laser(kat,'l1','n1',1)
space(kat,'s1','n1','n2',1)
mirror(kat,'m1','n2','n3',0.8,0.2)
space(kat,'s2','n3','n4',1)
mirror(kat,'m2','n4','n5',0.7,0.3)

photodiode(kat,'pd_cav','n4')
photodiode(kat,'pd_ref','n2')
photodiode(kat,'pd_trs','n5')

xaxis(kat, Scale.linear, [0,1e-6], kat.m2, kat.m2.xbeta, 100)

kat.maxtem = 1

[x,y_0,hdr] = kat.run(printout=0,printerr=0)

plot1D(x,y,hdr)
