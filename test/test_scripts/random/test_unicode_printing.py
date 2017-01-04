# -*- coding: utf-8 -*-
import pykat
from pykat import finesse

result = [] # prepare for the result
kat = finesse.kat() # create a fresh cat object
kat.verbose = False
kat.loadKatFile("LHO_IFO_maxtem2.kat") # load the conf

##############################################
## set excitationas, sensors, and some settings
##############################################
kat.parseCommands( "fsig sig1 ETMXHR 10 180")
kat.parseCommands( "fsig sig1 ETMYHR 10 0")
kat.parseCommands( "pd1 myomc 10  nOMC_HROC_trans")
kat.parseCommands( "xaxis sig1 f log 10 1k 10")
kat.parseCommands( "put myomc f1 $x1") # to follow
kat.parseCommands( "yaxis abs:deg")
kat.verbose = True
out = kat.run()    # do the computation
result.append(out['myomc'])   # append the result

print("PASSED")