__version__ = "0.6.2"

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

import finesse
import components
import detectors
import commands

from pykat.optics.gaussian_beams import beam_param




