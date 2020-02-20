import pykat
import numpy as np
from pykat.ifo import aligo
import scipy.constants as scc

QE   = 0.98
resp = 0.856 # A/W

base = aligo.make_kat()
base.IFO.remove_IMC_HAM2(True, False)

base.maxtem = 2
base.L0.P = 32

base.SRM.T   = 0.32
base.SR2.L    = 0.05
base.ETMX.L   = 35e-6
base.ETMY.L   = 35e-6

#https://alog.ligo-wa.caltech.edu/aLOG/index.php?callRep=47113
base.mod1.midx = 0.16
base.mod2.midx = 0.18

base.ITMY.setRTL(None, 1.42/100., 40e-6)
base.ITMX.setRTL(None, 1.5/100., 40e-6)
base.PRM.L.value /= 4

base.ITMY_lens.p = 1/base.ITMY_lens.f.value + 0e-6
base.ITMX_lens.p = 1/base.ITMX_lens.f.value + 0e-6

base.IFO.CARM.port.phase = -87 #92.5
base.IFO.PRCL.port.phase = -60 #119.8
base.IFO.PRCL.quad = 'I'

base.IFO.SRCL.port = base.IFO.POP_f2
base.IFO.SRCL.quad = 'I'

base.IFO.MICH.port.phase = 103.4
base.IFO.MICH.quad = 'Q'

mmx, mmy, _ = aligo.mismatch_cavities(base, 'ni')

base = aligo.setup(base, verbose=True)

base.DARM_lock.offset = -20e-3/(QE*resp) # Get 20mA of DC signal
base.SRCL_lock.gain /= 2
base.PRCL_lock.gain /= 2
base.DARM_lock.gain /= 2
base.CARM_lock.gain /= 2
base.MICH_lock.gain /= 2

base.SRCL_lock.accuracy /= 5
base.PRCL_lock.accuracy /= 5
base.DARM_lock.accuracy /= 5
base.CARM_lock.accuracy /= 5
base.MICH_lock.accuracy /= 5

base.IFO.zero_locks()

base.noxaxis = True

base.DARM_lock.offset = -20e-3/(QE*resp) # Get 20mA of DC signal
base.DARM_lock.accuracy /= 100

base.IFO.zero_locks()

AS_WFS = """
s s_OM3_ASL1 605m nOM3c nASL1a
lens ASL1 333m nASL1a nASL1b
s s_ASL1BS 0 nASL1b nASBSa
bs ASBS 0.5 0.5 0 0 nASBSa nASBSb nASBSc dump
s s_ASBS_AS_A 191m nASBSb nAS_A
s s_ASBS_AS_B 475m nASBSc nAS_B
"""

base.parse(AS_WFS)
base.parse(base.IFO.DHARD_P.fsig())
base.parse("pd2 AS_A_f2_Q_P_TF 45497355.0 110.0 $fs nAS_A")
base.parse("pdtype AS_A_f2_Q_P_TF y-split")