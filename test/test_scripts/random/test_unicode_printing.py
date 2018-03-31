# -*- coding: utf-8 -*-
import pykat
from pykat import finesse

result = [] # prepare for the result
kat = finesse.kat() # create a fresh cat object
kat.verbose = False
kat.load("LHO_IFO_maxtem2.kat") # load the conf

##############################################
## set excitationas, sensors, and some settings
##############################################
kat.parse( "fsig sig1 ETMXHR 10 180")
kat.parse( "fsig sig1 ETMYHR 10 0")
kat.parse( "pd1 myomc 10  nOMC_HROC_trans")
kat.parse( "xaxis sig1 f log 10 1k 10")
kat.parse( "put myomc f1 $x1") # to follow
kat.parse( "yaxis abs:deg")
kat.verbose = True
out = kat.run()    # do the computation
result.append(out['myomc'])   # append the result

print("PASSED")