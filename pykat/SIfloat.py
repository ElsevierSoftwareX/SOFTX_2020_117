import os
import re

#staticmethod
def SIfloat(value):
    if type(value)==list:
        return [convertToFloat(s) for s in value]
    else:
        return convertToFloat(value)
    
def convertToFloat(value):
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
    
    value = str(value)
    
    for i, j in __suffix.iteritems():
        value=value.replace(i, str(j))
        
    return float(value)
