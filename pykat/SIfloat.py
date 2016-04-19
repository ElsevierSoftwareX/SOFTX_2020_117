from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import re
import pykat.exceptions as pkex
import numpy as np

__suffix = {'y': 'e-24',  # yocto
            'z': 'e-21',  # zepto
            'a': 'e-18',  # atto
            'f': 'e-15',  # femto
            'p': 'e-12',  # pico
            'n': 'e-9',   # nano
            'u': 'e-6',   # micro
            'm': 'e-3',   # mili
            'c': 'e-2',   # centi
            'd': 'e-1',   # deci
            'k': 'e3',    # kilo
            'M': 'e6',    # mega
            'G': 'e9',    # giga
            'T': 'e12',   # tera
            'P': 'e15'   # peta
            }
            
__exceptions = ["$fs", "$mfs"]

def SIfloat(value):
    if value is None: 
        return value
    
    if str(value) in __exceptions:
        return value
        
    value = np.array(value)
    
    v = np.vectorize(convertToFloat)
    
    if value.size == 1:
        return float(v(value))
        
    a = v(value)
    
    if len(a) == 1:
        return a[0]
    else:
        return a
    
def convertToFloat(value):
    
    try:
        # first just try and convert the value
        return float(value)
        
    except ValueError as ex:
        # Catch any casting exeception
        value = value.strip()
        
        # only the last value can be an SI scaling letter 
        last = value[-1]

        if last in __suffix:
            # remove last character and append the SI scaling
            value = value[0:-1] + __suffix[last]
        else:
            raise pkex.BasePyKatException("Could not convert SI scaling in '{0}' to a float".format(value))
        
        try:   
            return float(value)
        except ValueError as ex:
            raise pkex.BasePyKatException("Unable to convert '{0}' into a float".format(value))
