from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

try:
    try:
        from ._version import __version__
    except ModuleNotFoundError as ex:
        __version__ = "develop"
except (NameError, ImportError) as ex:
    __version__ = "develop"

__min_req_finesse__ = '2.3.0'

# This flag is used to switch on the gui features in pkat at import time
USE_GUI = False
HAS_OPTIVIS = False

import six

########################
# Global helper functions
isContainer = lambda c: (not isinstance(c, six.string_types)) and hasattr(c, "__iter__")
is_container = lambda c: (not isinstance(c, six.string_types)) and hasattr(c, "__iter__")
########################

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
import pykat.style as style

from pykat.optics.gaussian_beams import BeamParam

from pykat.plotting import init_pykat_plotting
from pykat.style import use as set_plot_style

from .SIfloat import SIfloat

msg = "Could not find the finesse executable 'kat'" \
	"or you do not have the permissions to run it."

class nokat(object):
	def __getattribute__(self, attr):
		warn(msg)

try:
	kat = finesse.kat()
	v = kat.finesse_version()
except pkex.MissingFinesse:
	from warnings import warn
	warn(msg)
	kat = nokat()
	v = str(__min_req_finesse__)

version_found    = [int(_) for _ in v.split('-')[0].split('.')]
version_required = [int(_) for _ in __min_req_finesse__.split('-')[0].split('.')]

fulfilled = all([_ == __ for _,__ in zip(version_found, version_required)]) or any([_ > __ for _,__ in zip(version_found, version_required)])
            
if not fulfilled:
    raise pkex.BasePyKatException("Pykat %s requires Finesse version %s or higher. You have have %s" % (__version__ ,
                                                                                              str(__min_req_finesse__),
                                                                                              v))
def info():
    print("Pykat version: " + __version__)
    print("Pykat loaded from: " + __file__)
    if kat != nokat():
        print("Finesse version: " + str(v))
        print("Finesse loaded from: " + str(kat._finesse_exec()))
    else:
        print("Finesse could not be initialised")


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
