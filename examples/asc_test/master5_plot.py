from pykat import finesse
from pykat.commands import *
import pylab as pl
import numpy as np
import shelve
import sys

import matplotlib
formatter = matplotlib.ticker.EngFormatter(unit='', places=0)
formatter.ENG_PREFIXES[-6] = 'u'

import matplotlib.backends.backend_pdf
def printPDF(self, filename):
        pdfp = matplotlib.backends.backend_pdf.PdfPages('filename')
        pdfp.savefig(self,dpi=300,bbox_inches='tight')
        pdfp.close()

def main():
    print """
    --------------------------------------------------------------
    Example file for using PyKat to automate Finesse simulations
    Finesse: http://www.gwoptics.org/finesse
    PyKat:   http://www.gwoptics.org/pykat
    
    The file runs through the various pykat files which are used
    to generate the Finesse results reported in the document:
    `Comparing Finesse simulations, analytical solutions and OSCAR 
    simulations of Fabry-Perot alignment signals', LIGO-T1300345
    
    Run this file to plot the data generated with master5.py.
        
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
        
    print "--------------------------------------------------------"
    print " 12. Plotting beam tilt with thermal lens "
    gravity_tilt()

def gravity_tilt():
    xscale = 1.0
    yscale = 1.0e9
    data=np.loadtxt("thermal_gravity.txt")    
    maxtems = data[:,0]
    rows = len(maxtems)
    cols = len(data[0,:])
    fig=pl.figure()
    color_cycle=['b', 'c', 'r', 'k']
    lw=np.ones(rows+1)*3
    lw[-2]=2
    lw[-1]=1
    global data 
    global cols
    global lw
    for i in range(cols-1):
        #print i
        #pl.plot(maxtems,yscale*data[:,i+1],'-', color= color_cycle[i], linewidth=lw[i], label='ETM1')
        pl.scatter(maxtems,yscale*data[:,i+1],s=80,facecolors='none', edgecolors=color_cycle[i], label='xxx')


    pl.xlabel("maxtem")
    pl.ylabel("beam tilt [urad]")
    #pl.annotate('WFS1',xy=[0.42,70])
    #pl.annotate('WFS2',xy=[0.62,5])
    #pl.xlim([0,1])
    #pl.legend(loc=2)
    pl.grid()
    pl.draw()
    pl.show(block=0)
    printPDF(fig,'gravity_lens.pdf')
    

    
if __name__ == '__main__':
    main()
