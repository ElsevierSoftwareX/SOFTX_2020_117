from pykat import finesse
from pykat.commands import *
import copy
import pylab as pl
from collections import namedtuple
from collections import OrderedDict

import matplotlib
formatter = matplotlib.ticker.EngFormatter(unit='', places=0)
formatter.ENG_PREFIXES[-6] = 'u'

import matplotlib.backends.backend_pdf
def printPDF(self, filename):
	pdfp = matplotlib.backends.backend_pdf.PdfPages(filename)
	pdfp.savefig(self,dpi=300,bbox_inches='tight')
	pdfp.close()
	
def main():
	print """
	--------------------------------------------------------------
	Example file for using PyKat to automate Finesse simulations
	Finesse: http://www.gwoptics.org/finesse
	PyKat:   http://www.gwoptics.org/pykat
	
	The file runs through some Finesse simulations for the
	Glasgow speedmeter experiment, using the file:
	sagnac_base.kat
	
	Andreas Freise 30.10.2014
	--------------------------------------------------------------
	"""    
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
	extra.contents.append('qhdS sens 180 nout1 nout2')
	extra.contents.append('scale meter sens')
	
	kat.maxtem='off'
	Lambda=1064.0e-9
	result=OrderedDict()
	legend = {}
	
	# getting mass of the light mirror of cavity a
	global m
	
	m = kat.M1a.mass.value
	global f
	f = namedtuple('f', ('start','stop','points','data'))
	f.start=100
	f.stop = 2e4
	f.points = 100
	kat.parseKatCode('xaxis sig1 f log {0} {1} {2}'.format(f.start, f.stop, f.points-1))

	print "--------------------------------------------------------"
	print " 0. Run default file"
	out = kat.run()
	#f.data = np.logspace(np.log10(f.start), np.log10(f.stop), f.points)
	f.data=out.x # getting frequency vector from Finesse instead
	
	print "--------------------------------------------------------"
	print " 1. Computing SQL"
	hbar=6.62606957E-34/(2.0 *np.pi)
	SQL_x= np.sqrt( 4 * hbar / ( m * f.data**2 * 4 * np.pi * np.pi ))
	legend['SQL']=mylegend('SQL','k')
	result['SQL']=SQL_x
	result['default']=out.y
	legend['default']=mylegend('Stefan D. example file','b')

	print "--------------------------------------------------------"
	print " 2. Open ETM ports, i.e. replacing dump nodes"
	kat.M2a.remove()
	kat.M3a.remove()
	kat.parseCommands('bs1 M2a $T_ETM 0 0 0 nM2aw nM2an nM2aT1 nM2aT2')
	kat.parseCommands('bs1 M3a $T_ETM 0 0 0 nM3aw nM3an nM3aT1 nM3aT2')
	kat.M2b.remove()
	kat.M3b.remove()
	kat.parseCommands('bs1 M2b $T_ETM 0 0 0 nM2bw nM2bn nM2bT1 nM2bT2')
	kat.parseCommands('bs1 M3b $T_ETM 0 0 0 nM3bw nM3bn nM3bT1 nM3bT2')
	out=kat.run()	
	result['open']=out.y
	legend['open']=mylegend('Opening M2/M3 ports','r')

	print "--------------------------------------------------------"
	print " 3. Imbalanced BS"
	result['bs']=imbalanced_bs(kat)
	legend['bs']=mylegend('Imbalanced BS 49:51','g')

	print "--------------------------------------------------------"
	print " 3. Mass asymmetry"
	result['mass']=mass(kat)
	legend['mass']=mylegend('Mass asymmetry 10%','c')

	print "--------------------------------------------------------"
	print " Plotting results"
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
	
def plot_results(f, result, legend):
	fig=pl.figure()
	N=len(result)
	for i, T in enumerate(result.keys()):
		data=result[str(T)]
		pl.plot(f.data, data,legend[str(T)].lt, color=legend[str(T)].color, label=legend[str(T)].text, lw=legend[str(T)].lw)
		#line.set_dashes([12, 4]) 

	ax=pl.gca()
	ax.set_yscale('log')
	ax.set_xscale('log')
	pl.xlabel("f [Hz}")
	pl.ylabel("Sensitivity [m/sqrt(hz)]")
	pl.xlim([f.start,f.stop])
	pl.ylim([1E-20, 2E-17])
	pl.grid()
	pl.legend(loc=1)
	pl.draw()
	pl.show(block=0)
	printPDF(fig, "speedmeter_plots.pdf")
	
	
if __name__ == '__main__':
    main()

