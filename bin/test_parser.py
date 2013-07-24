import sys
sys.path.append("../")

from pykat import finesse
import pykat.profiling
import numpy as np
import pylab as pl

kat = finesse.kat()
#kat.load("D:\\finesse\\test\\kat_test\\physics\\mirror_astig_tilt_all_BH.kat")
kat.load('parse.kat')
kat.getPerformanceData = True

[run, perfdata] = kat.run(printout=0,printerr=0)

pykat.profiling.plotReducedPerformanceData(perfdata)


