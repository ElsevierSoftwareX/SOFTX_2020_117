from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from pykat import finesse
from pykat.commands import *
import pylab as pl
import numpy as np
import sys

import matplotlib
formatter = matplotlib.ticker.EngFormatter(unit='', places=0)
formatter.ENG_PREFIXES[-6] = 'u'

import matplotlib.backends.backend_pdf
def printPDF(self, filename):
        pdfp = matplotlib.backends.backend_pdf.PdfPages(filename)
        pdfp.savefig(self,dpi=300,bbox_inches='tight')
        pdfp.close()

def main():
    print("""
    --------------------------------------------------------------
    Example file for using PyKat to automate Finesse simulations
    Finesse: http://www.gwoptics.org/finesse
    PyKat:   http://www.gwoptics.org/pykat
        
    Run this file to plot the data generated with master4.py.
        
    Andreas Freise 16.01.2014
    --------------------------------------------------------------
    """)
    
    # shall we clear the workspace?
    # %reset -f
    # maybe close all plot windows?
    # close('all')
            
    print("--------------------------------------------------------")
    print(" Plotting beam tilt with thermal lens ")
    gravity_tilt()

    print("--------------------------------------------------------")
    print(" Plotting WFS signal with thermal lens ")
    asc_signal('asc_signals_5.txt', (0.3,0.15))
    asc_signal('asc_signals_50.txt', (0.3,0.15))

def asc_signal(filename,loc):
    xscale = 1.0
    yscale = 1.0
    global data 
    global cols
    global lw
    data=np.loadtxt(filename)
    # extracting only nonzero rows
    data = data[~np.all(data == 0, axis=1)]
    maxtems = data[:,0]
    [rows,cols] = data.shape
    fig=pl.figure()
    color_cycle=['b', 'c', 'r', 'k']
    labels=['ITM WFS1', 'ETM WFS1','ITM WFS2', 'ETM WFS2']
    lw=np.ones(rows+1)*3
    lw[-2]=2
    lw[-1]=1
    for i in [0,2,1,3]:
        pl.scatter(data[:,0],yscale*data[:,i+1],s=80,facecolors='none', edgecolors=color_cycle[i], label=labels[i]+"={0:.4} W/rad".format(data[-1,i+1]*yscale))
        pl.plot(data[:,0],yscale*data[:,i+1],'-',color=color_cycle[i], linewidth=1)

    pl.xlabel("maxtem")
    pl.ylabel("sensing matrix element [W/rad]")
    pl.xlim([0,max(maxtems)+1])
    #pl.ylim([-180,100])
    pl.legend(loc=loc)
    pl.grid()
    pl.draw()
    pl.show(block=0)
    pdfname=filename.split('.')[0]+'.pdf'
    printPDF(fig,pdfname)
    
    
def gravity_tilt():
    xscale = 1.0
    yscale = 1.0e9
    global data 
    global cols
    global lw
    data=np.loadtxt("thermal_gravity.txt")
    # extracting only nonzero rows
    data = data[~np.all(data == 0, axis=1)]
    maxtems = data[:,0]
    [rows,cols] = data.shape
    fig=pl.figure()
    color_cycle=['b', 'c', 'r', 'k']
    labels=['ITM WFS1', 'ETM WFS1','ITM WFS2', 'ETM WFS2']
    lw=np.ones(rows+1)*3
    lw[-2]=2
    lw[-1]=1
    for i in [0,2,1,3]:
        pl.scatter(data[:,0],yscale*data[:,i+1],s=80,facecolors='none', edgecolors=color_cycle[i], label=labels[i]+"={0:.4} nrad".format(data[-1,i+1]*yscale))
        pl.plot(data[:,0],yscale*data[:,i+1],'-',color=color_cycle[i], linewidth=1)

    pl.xlabel("maxtem")
    pl.ylabel("beam tilt [nrad]")
    pl.xlim([0,max(maxtems)+1])
    pl.ylim([-8,4])
    pl.legend(loc=3)
    pl.grid()
    pl.draw()
    pl.show(block=0)
    printPDF(fig,'gravity_lens.pdf')
    

    
if __name__ == '__main__':
    main()
