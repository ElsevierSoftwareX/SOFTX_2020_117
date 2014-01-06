from pykat import finesse
from pykat.commands import *
import pylab as pl
import numpy as np
import shelve
import copy
import sys
import shutil

from itertools import cycle
import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt
formatter = matplotlib.ticker.EngFormatter(unit='', places=0)
formatter.ENG_PREFIXES[-6] = 'u'


def main():

    print """
    --------------------------------------------------------------
    Example file for using PyKat to automate Finesse simulations
    Finesse: http://www.gwoptics.org/finesse
    PyKat:   https://pypi.python.org/pypi/PyKat/
    
    The file runs through the various pykat files which are used
    to generate the Finesse results reported in the document:
    `Comparing Finesse simulations, analytical solutions and OSCAR 
    simulations of Fabry-Perot alignment signals', LIGO-T1300345
    
    Run this file to create the data and master4_plot.py to plot 
    the results. Results are saved after each step and plots can
    be created at any time.
    
    Andreas Freise 06.12.2013
    --------------------------------------------------------------
    """
    
    # shall we clear the workspace?
    # %reset -f

    # making these global during testing and debugging
    #global kat
    #global out
    
    kat = finesse.kat(tempdir=".",tempname="test")
    kat.verbose = False
    
    tmpresultfile = 'myshelf2.dat'
    
    # loading data saved by master.py
    kat.loadKatFile('asc_base3.kat')
    try:
        tmpfile = shelve.open(tmpresultfile)
        result=tmpfile['result']
        tmpfile.close()
    except: raise Exception("Could not open temprary results file {0}".format(tmpresultfile))
        
    print "--------------------------------------------------------"
    print " 10. ASC signals for large misalignments (ETM)"
    asc_large(kat)
    

def asc_large(tmpkat):
    kat = copy.deepcopy(tmpkat)

    code_lock = """
    set err PDrefl_p re
    lock z $err 900 1p
    put* ETM phi $z
    noplot z
    """
        
    kat.parseKatCode(code_lock)
    kat.parseKatCode('yaxis abs')
    kat.parseKatCode('xaxis ETM ybeta lin 0 1u 100')
    maxtems = [1, 3, 7, 20]
    #kat.verbose=1
    xscale = 1e6
    yscale = 1e6
    global out
    tmpfilename = "datashelf2.dat"
    backupname = "datashelf2.dat.bck"
    out={}
    done_maxtems = []
    
    for tem in maxtems:
        done_maxtems.append(tem)
        print "  Calculating maxtem = %d " % tem
        kat.maxtem = tem
        out[str(tem)] = kat.run(printout=0,printerr=1)
        import os.path
        if os.path.isfile(tmpfilename):
            shutil.copyfile(tmpfilename, backupname)

        print " current results saved in: {0}".format(tmpfilename)
        tmpfile = shelve.open(tmpfilename)
        tmpfile['out']=out
        tmpfile['maxtems']=done_maxtems
        tmpfile.close()

def asc_signal(tmpkat):
    kat = copy.deepcopy(tmpkat)

    code_lock = """
    set err PDrefl_p re
    lock z $err 900 1p
    put* ETM phi $z
    noplot z
    """
    
    kat.parseKatCode(code_lock)
    kat.parseKatCode('yaxis abs')
    kat.noxaxis = True
    kat.maxtem=1

    signal=np.zeros((2, 2))
    kat.ITM.ybeta=1e-10
    kat.ETM.ybeta=0.0
    out = kat.run(printout=0,printerr=0)
    WFS1_idx=out.ylabels.index("WFS1_I")
    WFS2_idx=out.ylabels.index("WFS2_I")
    signal[0,0] = out.y[WFS1_idx]
    signal[1,0] = out.y[WFS2_idx]

    kat.ITM.ybeta=0.0
    kat.ETM.ybeta=-1e-10
    out = kat.run(printout=0,printerr=0)
    signal[0,1] = out.y[WFS1_idx]
    signal[1,1] = out.y[WFS2_idx]
    signal = signal *1e10
    sensors=('WFS1', 'WFS2')
    mirrors=('ITM', 'ETM')
    print "  ASC Matrix:"
    for i in range(2):
        print "  ", sensors[i], " ",
        for j in range(2):
            print "%12.10g" % signal[i,j],
        print mirrors[i]
    return signal
    
    
if __name__ == '__main__':
    main()
