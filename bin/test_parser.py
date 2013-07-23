import sys
sys.path.append("../")

from pykat import finesse
import numpy as np
import pylab as pl

kat = finesse.kat()
kat.load("parse.kat")

run = kat.run(printout=0,printerr=0)

pl.figure()

pl.plot(run.x,run.y)
pl.xlabel(run.xlabel)
pl.legend(run.ylabels)
pl.show()

