from pykat import finesse
from pykat.commands import *
import copy
from collections import namedtuple
from collections import OrderedDict

import matplotlib
BACKEND = 'Qt4Agg'
matplotlib.use(BACKEND)
import pylab as pl

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
	# getting the anle of incidence on end mirrors:
	AoI2 = float(kat.constants['AoI2'].value) / 180.*np.pi
    # fsig applied to a BS with AoI!=0 reduces the phase
	# change for a given signal in meters. Thus to make the
	# finesse results comaptible with the default SQL we need
	# to scale the results with cos(AoI):
	global AoIScale
	AoIScale=np.cos(AoI2)
	# setting frequency range
	global f
	f = namedtuple('f', ('start','stop','points','data'))
	f.start=100
	f.stop = 2e4
	f.points = 100
	kat.parseKatCode('xaxis sig1 f log {0} {1} {2}'.format(f.start, f.stop, f.points-1))

	# Reading Haixing Miao's reference data:
	datah1=np.loadtxt('QN_Sagnac_25ppm_loss.dat')
	datah2=np.loadtxt('QN_Sagnac_lossless.dat')
	
	print "--------------------------------------------------------"
	print " Run default file (no loss)"
	out = kat.run()
	#f.data = np.logspace(np.log10(f.start), np.log10(f.stop), f.points)
	f.data=out.x # getting frequency vector from Finesse instead
	
	print "--------------------------------------------------------"
	print " Computing SQL"
	hbar=6.62606957E-34/(2.0 *np.pi)
	SQL_x= np.sqrt( 4 * hbar / ( m * f.data**2 * 4 * np.pi * np.pi ))
	legend['SQL']=mylegend('SQL','k')
	result['SQL']=SQL_x
	result['default']=out.y*AoIScale
	legend['default']=mylegend('Finesse, no loss','m')

	result['H1']=datah1[:,1]
	legend['H1']=mylegend('Haixing, no loss','k')
	legend['H1'].lt='--.'

	print "--------------------------------------------------------"
	print " Run file with loss (or transmission) on end mirrors"
	L=12.5e-6
	#L=0
	#T=25e-6
	T=0
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
	result['loss']=out.y*AoIScale
	legend['loss']=mylegend('Finesse, 25ppm loss','b')

	result['H2']=datah2[:,1]
	legend['H2']=mylegend('Haixing, 25ppm loss','k')
	legend['H2'].lt='-.'
	
	print "--------------------------------------------------------"
	print " Additional imbalanced BS"
	result['bs']=imbalanced_bs(kat)*AoIScale
	legend['bs']=mylegend('Imbalanced BS 49:51','r')

	#print "--------------------------------------------------------"
	#print " Mass asymmetry"
	#result['mass']=mass(kat)*AoIScale
	#legend['mass']=mylegend('Mass asymmetry 10%','c')

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
	printPDF(fig, "speedmeter_plots.pdf")
	
	
if __name__ == '__main__':
    main()

