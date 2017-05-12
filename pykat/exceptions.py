from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pykat.external.six as six
if six.PY2:
	import exceptions
import os, sys

def PrintError(message, exception=None):
    size = 62
    
    print("\033[91m")
    
    try:
        from textwrap import wrap, fill
        print ("-" * size)
        
        for a in wrap(message, size):
            print(a)
            
        if exception is not None:
            for a in wrap(str(exception.msg), size):
                a = a.replace("*** ", "\n")
                a = a.replace("** ", "\n")
                print(a)
            
        print ("-" * size)
    finally:
        print ("\033[0m")
        sys.exit(1)

def printWarning(message, exception=None):
    size = 62
    
    print("\033[93m")
    
    try:
        from textwrap import wrap, fill
        
        for a in wrap(message, size):
            print(a)

        if exception is not None:
            if hasattr(exception, "msg"):
                msg = str(exception.msg)
            else:
                msg = str(exception)
                
            for a in wrap(msg, size):
                a = a.replace("*** ", "\n")
                a = a.replace("** ", "\n")
                print(a)
            
    finally:
        print ("\033[0m")
    

class BasePyKatException(Exception):
    def __init__(self, msg):
        self.msg = msg
        
    def __str__(self):
        return self.msg

class FinesseParse(BasePyKatException) :    
    def __init__(self, msg):
        BasePyKatException.__init__(self, "Error parsing Finesse input\n{0}".format(msg))
    
class MissingFinesseEnvVar(BasePyKatException) :    
    def __init__(self):
        BasePyKatException.__init__(self, "The environment variable FINESSE_DIR was not defined")

class MissingFinesse(BasePyKatException) :    
    def __init__(self):
        BasePyKatException.__init__(self, "Could not find the finesse executable 'kat' in '{0}'," \
                                     "or you do not have the permissions to run it." \
                                      .format(os.environ.get('FINESSE_DIR')))
    
class FinesseRunError(BasePyKatException) :
    def __init__(self, err, kat):
        self.__err = err
        self.__kat = kat
        
        BasePyKatException.__init__(self, self.__err)
        #BasePyKatException.__init__(self, "{0}".format(self.__err))
        
