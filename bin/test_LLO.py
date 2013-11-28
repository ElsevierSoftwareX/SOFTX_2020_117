from pykat import finesse
from pykat.commands import xaxis
import pylab as pl
import numpy as np
from pykat_LLO import LLO

kat = finesse.kat()

kat.parseCommands(LLO.Laser)
kat.parseCommands(LLO.PR)
kat.parseCommands(LLO.BS)
kat.parseCommands(LLO.X_ARM)

kat.HRBS.r_ap = 37e-2

kat.parseCommands("pd pdIN nin")
kat.noxaxis = True

out = kat.run()
