__version__ = "0.6.2"

# This flag is used to switch on the gui features in pkat at import time
USE_GUI = False

import pykat.exceptions as pkex

NoGUIException = pkex.BasePyKatException("No PyQt4 module was found so cannot open a GUI")

import pykat.finesse as finesse
import pykat.components as components
import pykat.detectors as detectors
import pykat.commands as commands

from pykat.optics.gaussian_beams import beam_param




