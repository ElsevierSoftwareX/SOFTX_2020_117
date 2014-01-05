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

import matplotlib.backends.backend_pdf
def printPDF(self):
        filename = 'test.pdf'
        pdfp = matplotlib.backends.backend_pdf.PdfPages('large_ITM.pdf')
        pdfp.savefig(self,dpi=300,bbox_inches='tight')
        pdfp.close()

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
    
    Run this file to plot the data generated with master3.py.
        
    Andreas Freise 06.12.2013
    --------------------------------------------------------------
    """
    
    # shall we clear the workspace?
    # %reset -f
    # maybe close all plot windows?
    # close('all')
    
    # making these global during testing and debugging
    #global kat
    #global out
    #global result
        
    tmpresultfile = 'myshelf2.dat'
    
    # loading data saved by master.py
    try:
        tmpfile = shelve.open(tmpresultfile)
        result=tmpfile['result']
        tmpfile.close()
    except: raise Exception("Could not open temprary results file {0}".format(tmpresultfile))
                
    print "--------------------------------------------------------"
    print " 9. ASC signals for large misalignments"
    asc_large()
    

def asc_large():
    xscale = 1e6
    yscale = 1e6
    tmpfilename = "datashelf1.dat"
    backupname = "datashelf1.dat.bck"

    try:
        tmpfile = shelve.open(tmpfilename)
        out=tmpfile['out']
        maxtems=tmpfile['maxtems']
        tmpfile.close()
    except: raise Exception("Could not open temprary results file {0}".format(tmpfilename))

    fig=pl.figure()
    color_cycle=['b', 'c', 'r', 'k']
    N=len(maxtems)
    lw=np.ones(N)*3
    lw[-2]=2
    lw[-1]=1
    for i, tem in zip(range(len(maxtems)), maxtems):
        data=out[str(tem)]
        WFS1_idx=data.ylabels.index("WFS1_I")
        WFS2_idx=data.ylabels.index("WFS2_I")
        pl.plot(xscale*data.x,yscale*data.y[:,WFS1_idx],'-', color= color_cycle[i], linewidth=lw[i], label='maxtem {0}'.format(tem))
        line, = pl.plot(xscale*data.x,yscale*data.y[:,WFS2_idx],'-', color = color_cycle[i], linewidth=lw[i])
        #line.set_dashes([12, 4]) 

    osc1=np.loadtxt("OSCAR_large_tilt_ITM.txt",comments='%')

    x=xscale*osc1[:,0]
    y=yscale*osc1[:,1]
    pl.scatter(x,y,s=80,facecolors='none', edgecolors='k', label='OSCAR')
    y=yscale*osc1[:,2]
    pl.scatter(x,y,s=80,facecolors='none', edgecolors='k')
    pl.xlabel("ITM ybeta [urad]")
    pl.ylabel("Alignment signal [uW]")
    pl.annotate('WFS1',xy=[0.42,70])
    pl.annotate('WFS2',xy=[0.62,5])
    pl.xlim([0,1])
    pl.legend(loc=2)
    pl.grid()
    pl.draw()
    pl.show(block=0)
    printPDF(fig)
    

    
if __name__ == '__main__':
    main()
