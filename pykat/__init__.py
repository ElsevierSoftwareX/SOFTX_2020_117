from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

__version__ = "1.0.18"

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

from .SIfloat import SIfloat


SI = {'yocto': 1E-24,  # yocto
    'zepto': 1E-21,  # zepto
    'atto': 1E-18,  # atto
    'femto': 1E-15,  # femto
    'pico': 1E-12,  # pico
    'nano': 1E-9,   # nano
    'micro': 1E-6,   # micro
    'milli': 1E-3,   # mili
    'centi': 1E-2,   # centi
    'deci': 1E-1,   # deci
    None: 1E-0,  
    'kilo': 1E3,    # kilo
    'mega': 1E6,    # mega
    'giga': 1E9,    # giga
    'tera': 1E12,   # tera
    'peta': 1E15   # peta
    }

SIlabel = {'yocto': 'y',# yocto
           'zepto': 'z',# zepto
            'atto': 'a',# atto
           'femto': 'f',# femto
            'pico': 'p',# pico
            'nano': 'n',# nano
            'micro':'u',# micro
            'milli':'m',# mili
            'centi':'c',# centi
            'deci': 'd',# deci
            None:   'k',
            'kilo': 'M',# kilo
            'mega': 'G',# mega
            'giga': 'T',# giga
            'tera': 'P',# tera
            'peta': 'y',# peta
    }