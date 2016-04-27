"""
This example is a new feature regarding the use of Finesse as a shared library rather than 
as an executable. It is still very experimental! The shared library is accessed using
multiprocessing, which spawns an entire new process to run Finesse in. To use this requires
firstly generating pylibkat.py, which involves generating a ctypes wrapper using 
https://code.google.com/p/ctypesgen/. The command I use to generate the wrapper is:

    /Users/ddb/svn/ctypesgen-read-only/ctypesgen.py -o pylibkat.py ./kat.h 

You then also need to build libkat.so shared library from the python_shared branch of Finesse.
Once the above is done the following script may or may not work...
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


import ctypes
import pylibkat
import scipy
from scipy.optimize import minimize
import pylab
import time
import numpy as np
from multiprocessing import Process, Manager, Value
import pykat

def callback(lkat, maxphi):
    """
    This callback will be run in a completly different process to that which this
    script started. So we can't simply pass values back and forth. The callback
    arguments are:
        lkat - The handle to the finesse instance
        maxphi - a shared object between this and the main process, see multiprocessing.Value
    """
    print ("Entering callback...")
    
    # first we need to get a handle on the internals of Finesse
    inter = pylibkat.interferometer.in_dll(lkat, "inter")
    
    m1 = inter.mirror_list[0]
    m2 = inter.mirror_list[1]   

    circ = inter.output_data_list[0]
    
    def Fmin(x):
        # change any variables we are interested in
        m1.phi = x[0]
        # now step the simulation forward
        lkat._pykat_step()
        # return our cost function
        return -1 * circ.re       

    minimize(Fmin, [maxphi.value], method='Nelder-Mead')

    maxphi.value = m1.phi
    
    print ("Process: Maximum power =", circ.re)
    print ("Process: Mirror tuning =", m1.phi)

cmd = """
l l1 1 0 n1
s s1 1 n1 n2
m m1 0.99 0.01 0 n2 n3
s s2 100 n3 n4
m m2 0.99 0.01 0 n4 n5
pd circ n3

noxaxis
maxtem 2

attr m1 Rc -1000
attr m2 Rc 1000 
cav c1 m1 n3 m2 n4
"""

kat = pykat.finesse.kat()

kat.parseCommands(cmd)

maxphi = Value('d', 0)

p = kat.getProcess(callback, maxphi=maxphi)
p.start()
p.join()

print ("Host:    Received maximum phi =", maxphi.value)
