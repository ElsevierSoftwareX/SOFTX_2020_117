import os
import re

#staticmethod
def SIfloat(value):
    if type(value)==list:
        return [convertToFloat(s) for s in value]
    else:
        return convertToFloat(value)
    
def convertToFloat(value):
    __prefix = {'y': 1e-24,  # yocto
                'z': 1e-21,  # zepto
                'a': 1e-18,  # atto
                'f': 1e-15,  # femto
                'p': 1e-12,  # pico
                'n': 1e-9,   # nano
                'u': 1e-6,   # micro
                'm': 1e-3,   # mili
                'c': 1e-2,   # centi
                'd': 1e-1,   # deci
                'k': 1e3,    # kilo
                'M': 1e6,    # mega
                'G': 1e9,    # giga
                'T': 1e12,   # tera
                'P': 1e15,   # peta
                'E': 1e18,   # exa
                'Z': 1e21,   # zetta
                'Y': 1e24,   # yotta
                }
    value = str(value)
    for i, j in __prefix.iteritems():
        value=value.replace(i, str(j))
    return float(value)
