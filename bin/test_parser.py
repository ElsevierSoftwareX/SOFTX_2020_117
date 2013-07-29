import sys
sys.path.append(".")

from pykat import finesse, profiling
import numpy as np
import pylab as pl

kat = finesse.kat()
#kat.load("D:\\finesse\\test\\kat_test\\physics\\mirror_astig_tilt_all_BH.kat")
kat.load('parse.kat')
kat.getPerformanceData = True

[run, perfdata] = kat.run(printout=0,printerr=0)

l,t,fig = profiling.plotReducedPerformanceData(perfdata)

fig.savefig("Timing.png",pad_inches=0.1, bbox_inches='tight')
fig.show()
    