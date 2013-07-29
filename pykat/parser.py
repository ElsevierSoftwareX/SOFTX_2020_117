import os
import exceptions
import numpy as np

from pykat.node_network import NodeNetwork
from pykat.detectors import Detector
from pykat.components import Component
from pykat.commands import Command, xaxis

components = np.array(['m','m1','m2','l','s','bs','bs1','bs2','pd','pd*'])
commands = np.array(['attr','tem','tem*','gauss','gauss*','gauss**','cav','conf'])


# some commands we ignore, we do the plotting with pyhton
# so don't need 
ignore = ['gnuterm']
def parse_kat_file(kat_filename):
    
    kat_cmps = [] # holds the components found in kat file
    kat_cmds = [] # holds the commands found in kat file
    
    katfile = open(kat_filename,'r')
    
    for line in katfile.readlines():
        arg = line.split(' ')[0]
        
        # c
        if (components == arg).any():
            kat_cmps.append()
            
        elif (commands == arg).any():
            print ""

def parse_m(line):
    return line
    