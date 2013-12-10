# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 14:18:17 2013

@author: Sean
"""

import sys
sys.path.append("../")

import pykat
from pykat.utilities.optics.gaussian_beams import gauss_param
import pykat.finesse as finesse
from pykat.commands import xaxis
import pylab as pl
import numpy as np
import math

code = """
%------------------------------------------------------------------------
% Finesse input file to plot the phase of light field reflected from a
% beam splitter to show the way lengths and positions are handled
% Andreas Freise 15.08.2009
%------------------------------------------------------------------------
                 
l l1 1 0 n1    % laser with P=1W at the default frequency
s s1 1 1 n1 n2 % space of 1m length
bs b1 1 0 0 0 n2 n3 dump dump % beam splitter as `turning mirror'
s s2 1 1 n3 n4 % another space of 1m length
ad ad1 0 n4     % amplitude detector
 
% for the plot we perform two sequenctial runs of Finesse
% 1) first trace: change microscopic position of beamsplitter
xaxis b1 phi lin 0 180 100
% 2) second trace: change length of space s1
% xaxis s1 L lin 1 2 100
 
yaxis deg     % plotting the phase of the results
"""

kat = finesse.kat()
kat.parseCommands(code)

maxtem = np.arange(0, 2, 2)

for tem in maxtem:
    print "Calculating maxtem ", tem, "..."
    kat.maxtem = tem
    r = kat.run()
    pl.plot(r.x, r.y, label="maxtem={0}".format(tem))

pl.ylabel("Phase [deg]")
pl.xlabel("Tuning [deg]")
pl.legend()
pl.show()
    
