import matplotlib
matplotlib.use('Agg')
import pykat
import json
#import matplotlib.pyplot as plt

import pykat.ifo.aligo as aligo
import pykat.ifo.aligo.plot

kat = aligo.make_kat()


kat.L0.P.value = 1
kat.verbose = False

# We'll be regenerating these blocks so remove them now
kat.removeBlock('locks', failOnBlockNotFound=False)
kat.removeBlock('errsigs', failOnBlockNotFound=False)
kat.removeBlock('powers', failOnBlockNotFound=False)
    
kat.ETMX.mass = None
kat.ETMY.mass = None
kat.ITMX.mass = None
kat.ITMY.mass = None

kat.maxtem = 0
kat.phase  = 2

kat.IFO.adjust_PRC_length()

###############################################################
# Pretuning ensures each of the ARMs, PRC and SRC are on the correct resonance to being with
# Run a copy of aLIGO without the modulators to get the pretuning
_ = kat.deepcopy()
_.IFO.remove_modulators()
aligo.pretune(_, verbose=True)

# Apply the tunings to our base kat file
kat.IFO.apply_tunings(_.IFO.get_tunings())

aligo.pretune_status(kat)
kat.IFO.lengths_status()
###############################################################


###############################################################
# Set DC offset and plot
DCoffset = 20e-12 / kat.lambda0 * 180.0 # 20 pm
kat.IFO.set_DC_offset(DCoffset=DCoffset, verbose=True)

#aligo.plot.error_signals(kat, xlimits=[-1e-2, 1e-2])
###############################################################


###############################################################
errsigs_cmds = kat.IFO.add_errsigs_block()

# Generates a dictionary of the lock values to use
locks = aligo.generate_locks(kat, verbose=True)

# Takes these values and then generates the commands and adds them to
# the lock block of the kat
lock_cmds = kat.IFO.add_locks_block(locks, verbose=True)

# Running locks once
_ = kat.deepcopy()
_.noxaxis = True
_.verbose = True
out = _.run()

kat.IFO.apply_lock_feedback(out)
###############################################################

print("\33[92mDONE\33[0m")
