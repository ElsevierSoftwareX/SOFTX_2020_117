from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

__version__ = "1.0.19"

# This flag is used to switch on the gui features in pkat at import time
USE_GUI = False
HAS_OPTIVIS = False

import imp

try:
	imp.find_module('optivis')
	HAS_OPTIVIS = True
except ImportError:
	HAS_OPTIVIS = False

import pykat.exceptions as pkex

NoGUIException = pkex.BasePyKatException("No PyQt4 module was found so cannot open a GUI")
    
import pykat.finesse as finesse
import pykat.components as components
import pykat.detectors as detectors
import pykat.commands as commands

from pykat.optics.gaussian_beams import BeamParam

from pykat.plotting import init_pykat_plotting

