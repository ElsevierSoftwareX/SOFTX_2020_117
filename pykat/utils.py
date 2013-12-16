# -*- coding: utf-8 -*-
"""
Created on Sat Feb 02 12:17:50 2013

@author: Daniel
"""
from pykat.components import *

def isSpace(obj):
    if obj == None:
        return False
    else:
        return isinstance(obj,space)