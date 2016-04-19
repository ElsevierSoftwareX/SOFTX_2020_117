from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from pykat import finesse
from pykat.commands import *
import pylab as pl
import numpy as np
import pickle
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
	PyKat:	 http://www.gwoptics.org/pykat
		
	Run this file to plot the data generated with master3.py.
		
	Andreas Freise 16.01.2014
	--------------------------------------------------------------
	""")
	
	# shall we clear the workspace?
	# %reset -f
	# maybe close all plot windows?
	# close('all')
	
	# making these global during testing and debugging
	#global kat, out, result
		
	print("--------------------------------------------------------")
	print(" Plotting ASC signals for large misalignments")
	asc_large('ITM')
	asc_large('ETM')
	

def asc_large(mir_name):
	xscale = 1e6
	yscale = 1e6
	if mir_name == 'ITM':
		ysign = 1.0
	else:
		ysign = -1.0
		
	tmpfilename = "datashelf_{0}.dat".format(mir_name)
	backupname = "datashelf_{0}.dat.bck".format(mir_name)

	try:
		readpickle = pickle.load(open(tmpfilename, "rb"))
		out=readpickle['out']
		maxtems=readpickle['maxtems']
	except: raise Exception("Could not open temprary results file {0}".format(tmpfilename))

	fig=pl.figure()
	color_cycle=['b', 'c', 'r', 'k']
	N=len(maxtems)
	lw=np.ones(N)*3
	lw[-2]=2
	lw[-1]=1
	for i, tem in enumerate(maxtems):
		data=out[str(tem)]
		WFS1_idx=data.ylabels.index("WFS1_I")
		WFS2_idx=data.ylabels.index("WFS2_I")
		pl.plot(xscale*data.x,ysign*yscale*data.y[:,WFS1_idx],'-', color= color_cycle[i], linewidth=lw[i], label='maxtem {0}'.format(tem))
		line, = pl.plot(xscale*data.x,ysign*yscale*data.y[:,WFS2_idx],'-', color = color_cycle[i], linewidth=lw[i])
		#line.set_dashes([12, 4]) 

	osc1=np.loadtxt("OSCAR_large_tilt_{0}.txt".format(mir_name),comments='%')

	x=xscale*osc1[:,0]
	y=yscale*osc1[:,1]
	pl.scatter(x,y,s=80,facecolors='none', edgecolors='k', label='OSCAR')
	y=yscale*osc1[:,2]
	pl.scatter(x,y,s=80,facecolors='none', edgecolors='k')
	pl.xlabel("{0} vertical misalignment [urad]".format(mir_name))
	pl.ylabel("Alignment signal [uW]")
	if mir_name == 'ITM':
		pl.annotate('WFS1',xy=[0.42,70])
		pl.annotate('WFS2',xy=[0.62,5])
		pl.legend(loc=2)
	else:
		pl.annotate('WFS1',xy=[0.42,10])
		pl.annotate('WFS2',xy=[0.62,-70])
		pl.legend(loc=3)
	pl.xlim([0,1])
	pl.grid()
	pl.draw()
	pl.show(block=0)
	printPDF(fig, "large_{0}.pdf".format(mir_name))
	

if __name__ == '__main__':
	main()
