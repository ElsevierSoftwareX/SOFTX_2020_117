# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 09:56:53 2013

PyKat - Python interface and wrapper for FINESSE
Copyright (C) 2013 Daniel David Brown

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

Contact at ddb@star.sr.bham.ac.uk

@author: Daniel Brown
"""

from IPython.parallel import Client
import sys
import os

def _run(commands, pwd):
    import os
    os.chdir(pwd)
    
    import pykat

    kat = pykat.finesse.kat()
    kat.parseCommands(commands)
    out = kat.run(rethrowExceptions=True)
    
    return out

class parakat(object):
    """
    Uses the ipython clustering for running kat objects in parallel.
    
    To use this you must have started an ipython cluster on your computer.
    From a new terminal use the command:
        
        ipcluster start -n 4
        
    This will start a cluster with 4 workers.
    
    To run a kat object use:
    
        pk = parakat()
        pk.run(kat1)
        pk.run(kat2)
        pk.run(kat3)
        
        outs = pk.getResults()
    
    The list 'outs' will contain the katRun object you'd normal get if you 
    had just called, kat1.run(), etc. The results list is matched to order
    in which you run the kats.
    
    If you need to stop long running kat processes the chances are you will
    also need to kill the ipython cluster process, as sometimes they carry
    on running.
    """
    
    def __init__(self):
        self._rc = Client()
        self._lview = self._rc.load_balanced_view()
        self._lview.block = False
        self._results = []
        
    def run(self, kat):
        print(kat)
        self._results.append(self._lview.apply_async(_run, "".join(kat.generateKatScript()), os.getcwd()))
        print(self._results)
        
    def getResults(self):
        out = []
        
        self._lview.wait(self._results)
        
        for done in self._results:
            out.append(done.get())
            
        return out
