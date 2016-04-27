from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from pykat import finesse
from pykat.commands import *
import copy
import pylab as pl
from collections import namedtuple
from collections import OrderedDict

import matplotlib
formatter = matplotlib.ticker.EngFormatter(unit='', places=0)
formatter.ENG_PREFIXES[-6] = 'u'

global fig_sagnac

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
	
	The file runs through some Finesse simulations for the
	Glasgow speedmeter experiment, using the file:
	sagnac_base.kat
	
	Andreas Freise 30.10.2014
	--------------------------------------------------------------
	""")    
    # defining variables as global for debugging
	global kat
	global out
	global result
	global legend
    
	# for debugging we might need to see the temporay file:
	kat = finesse.kat(tempdir=".",tempname="test")
	kat.verbose = False
	
	kat.loadKatFile('sagnac_base.kat')
	extra = kat._kat__blocks['NO_BLOCK']

	# adding homodyne detector
	#kat.addLine('qhdS sens 180 nout1 nout2')
	extra.contents.append('qhdS sens 180 nout1 nout2')
	extra.contents.append('scale meter sens')
	
	kat.maxtem='off'
	Lambda=1064.0e-9
	result=OrderedDict()
	legend = {}
	
	# getting mass of the light mirror of cavity a
	global m
	
	m = kat.M1a.mass.value
	AoI2 = float(kat.constants['AoI2'].value) / 180.*np.pi
	global f
	f = namedtuple('f', ('start','stop','points','data'))
	f.start=100
	f.stop = 2e4
	f.points = 100
	kat.parseKatCode('xaxis sig1 f log {0} {1} {2}'.format(f.start, f.stop, f.points-1))

	# Reading Haixing Miao's reference data:
	datah1=np.loadtxt('QN_Sagnac_lossless.dat')
	datah2=np.loadtxt('QN_Sagnac_25ppm_loss.dat')

	# Reading Stefan D. example data:
	data=np.loadtxt('Stefan_data.txt')
	
	print ("--------------------------------------------------------")
	print (" Run default file (no loss)")
	out = kat.run()
	#f.data = np.logspace(np.log10(f.start), np.log10(f.stop), f.points)
	f.data=out.x # getting frequency vector from Finesse instead
	
	print ("--------------------------------------------------------")
	print (" Computing SQL")
	hbar=6.62606957E-34/(2.0 *np.pi)
	SQL_x= np.sqrt( 4 * hbar / ( m * f.data**2 * 4 * np.pi * np.pi ))
	legend['SQL']=mylegend('SQL','k')
	result['SQL']=SQL_x

	result['H1']=datah1[:,1]
	legend['H1']=mylegend('Haixing, no loss','g')
	legend['H1'].lt='--.'
	legend['H1'].lw=2

	result['default']=out.y*np.cos(AoI2)
	legend['default']=mylegend('no loss','m')


	
	print ("--------------------------------------------------------")
	L=0
	T=12.5e-6
	R=1-T-L
	kat.M2a.R=R
	kat.M3a.R=R
	kat.M2b.R=R
	kat.M3b.R=R
	kat.M2a.T=T
	kat.M3a.T=T
	kat.M2b.T=T
	kat.M3b.T=T
	kat.M2a.L=L
	kat.M3a.L=L
	kat.M2b.L=L
	kat.M3b.L=L
	out=kat.run()

	result['H2']=datah2[:,1]
	legend['H2']=mylegend('Haixing, 25ppm loss','g')
	legend['H2'].lt='-.'
	legend['H2'].lw=5
	
	result['loss']=out.y*np.cos(AoI2)
	legend['loss']=mylegend('25ppm loss','b')


	
	#result['S_sym']=data[:,1]
	#legend['S_sym']=mylegend('Stefan D, balanced','k')
	#legend['S_sym'].lt='-.'

	print ("--------------------------------------------------------")
	print (" 3. Imbalanced BS")
	#result['bs']=imbalanced_bs(kat)*np.cos(AoI2)
	#legend['bs']=mylegend('Imbalanced BS 49:51','r')

	#result['S_imb']=data[:,2]
	#legend['S_imb']=mylegend('Stefan D, imbalanced','k')
	#legend['S_imb'].lt='--.'

	print ("--------------------------------------------------------")
	print (" 3. Mass asymmetry")
	#result['mass']=mass(kat)
	#legend['mass']=mylegend('Mass asymmetry 10%','c')

	print ("--------------------------------------------------------")
	print (" Plotting results")
	plot_results(f, result, legend)


# custom class to store text, color, line type and line width for
# legend entries
class mylegend():
	def __init__(self, _text, _color):
		self.text=_text
		self.color=_color;
	lw=2
	lt='-'
	
def mass(tmpkat):
	kat = copy.deepcopy(tmpkat)
	kat.M1a.mass=m*0.9
	kat.M1b.mass=m*1.1
	out=kat.run()
	return out.y
	
	
def imbalanced_bs(tmpkat):
	kat = copy.deepcopy(tmpkat)
	kat.M6.R=0.49
	kat.M6.T=0.51
	out=kat.run()
	return out.y


def update_plot(f, result):
	for i, T in enumerate(result.keys()):		
		data=result[str(T)]
		for j,f1 in enumerate(f):
			lines[i].set_ydata(data(j))
			fig.canvas.draw()
		
def plot_results(f, result, legend):
	fig=pl.figure(232)
	fig.clear()
	N=len(result)
	lines={}
	for i, T in enumerate(result.keys()):
		data=result[str(T)]
		lines[i],=pl.plot(f.data, data,legend[str(T)].lt, color=legend[str(T)].color, label=legend[str(T)].text, lw=legend[str(T)].lw)
		#line.set_dashes([12, 4]) 

	ax=pl.gca()
	ax.set_yscale('log')
	ax.set_xscale('log')
	pl.xlabel("f [Hz]")
	pl.ylabel("Sensitivity [m/sqrt(Hz)]")
	pl.xlim([f.start,f.stop])
	pl.ylim([1E-20, 2E-17])
	pl.grid()
	pl.legend(loc=1)
	pl.draw()
	pl.show(block=0)
	#printPDF(fig, "speedmeter_plots.pdf")
	
	
if __name__ == '__main__':
    main()

