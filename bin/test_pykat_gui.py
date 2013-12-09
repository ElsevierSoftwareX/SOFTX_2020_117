import sys
sys.path.append("../")

from pykat import finesse
from pykat.detectors import *
from pykat.components import *
from pykat.commands import *
from pykat.structs import *
import numpy as np
import pylab as pl

code = """
#l l1 1 0 0 n1

s s1 1 n3 n4
"""

kat = finesse.kat(kat_code = code)

kat.openGUI()


