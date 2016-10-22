from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

import copy
from collections import namedtuple
from collections import OrderedDict
import shelve

import pylab as pl

import pykat
from pykat.components import *

from pykat.external.progressbar import ProgressBar, ETA, Percentage, Bar, Timer
from pykat.optics.maps import *
from pykat.optics.gaussian_beams import HG_mode, beam_param
from pykat.optics.fft import *
#from pykat.tools.plotting.tools import plot_field, plot_propagation

from aligo import *

def main():
	print("""
	--------------------------------------------------------------
	Example file for using PyKat for an FFT-based simulation
	of an Advanced LIGO arm cavity
	PyKat:   http://www.gwoptics.org/pykat
	Advanced LIGO: http://www.advancedligo.mit.edu/

	Requires surface map file: etm08_virtual.txt
	Andreas Freise 20.12.2014
	--------------------------------------------------------------
	""") 
    # defining variables as global for debugging

	tmpresultfile = 'myshelf1.dat'
	result={}

	#Advanced LIGO parameters
	# length [m] 3994.5
	# FSR [kHz] 37.5
	# T1 1.4%
	# T2 5ppm
	# finesse 445
	# FWHM [Hz] 84
	# mirror diameter: 0.32 m
	# ITM RC [m] 1934
	# ETM RC [m] 2245
	# w0 [cm] 1.2
	# w1 [cm] 5.3
	# w2 [cm] 6.2
	# z1 [m] -1834

	# loading kat file to get parameters and to compute input beam parameters
	global kat, out
	kat = pykat.finesse.kat()
	kat.verbose = False
	kat.loadKatFile('aligo_Xarm.kat')

	# setting ITM T to larger value for better plots
	kat.itmX.T=0.1
	kat.itmX.R=0.9
	Lambda = kat.lambda0
	LX=kat.LX.L.value
	kat.maxtem=0
	out = kat.run()
	w0=out.y[0][0]
	z0=-out.y[0][1]
		
	# load and create mirror maps
	global itm, etm
	surface=read_map('etm08_virtual.txt')
	itm=curvedmap('itm_Rc',surface.size,surface.step_size, -1.0*abs(kat.itmX.Rc.value))
	etm=curvedmap('etm_Rc',surface.size,surface.step_size, -1.0*abs(kat.etmX.Rc.value))
	#itm.plot()
	#etm.plot()
	# apply measured map to etm, using 20 times larger distortions
	etm.data = etm.data + surface.data*surface.scaling/etm.scaling*20

	# setup grid for FFT propagation
	[xpoints,ypoints] = surface.size
	xsize = xpoints * surface.step_size[0]
	ysize = ypoints * surface.step_size[1]
	xoffset = 0.0
	yoffset = 0.0

	global shape
	shape = grid(xpoints, ypoints, xsize, ysize, xoffset, yoffset)
	x = shape.xaxis
	y = shape.yaxis
	result['shape']=shape

	# generate roughly mode-matched input beam
	global laser
	gx = beam_param(w0=w0, z=z0)
	beam = HG_mode(gx,gx,0,0)
	laser = beam.Unm(x,y) 

	# some debugging plots
	"""
	plot_field(laser)
	Lrange= np.linspace(0,4000,200)
	plot_propagation(laser, shape, Lambda, 0, 1, Lrange, 1)
	laser1=FFT_propagate(laser,shape,Lambda,LX,1)
	laser2=np.sqrt(kat.etmX.R.value)*FFT_apply_map(laser1, etm, Lambda)
	laser3=FFT_propagate(laser2,shape,Lambda,LX,1)
	Lrange= np.linspace(0,4000,200)
	plot_propagation(laser2, shape, Lambda, 0, 1, Lrange, 1)
	plot_field(laser3)
	"""

	precompute_roundtrips(shape, laser, kat)
	
	# now save any `result' variables:
	tmpfile = shelve.open(tmpresultfile)
	tmpfile['result']=result
	tmpfile.close()

		
def precompute_roundtrips(shape, laser, kat):
	Lambda=kat.lambda0
	LX=kat.LX.L.value
	R=kat.etmX.R.value*kat.itmX.R.value
	Loss = 1-R
	accuracy=100E-6
	print("cavity loss: {0}".format(Loss))	
	N=int(required_roundtrips(Loss,accuracy))
	print("required rountrips: {0} (for accuracy of {1})".format(N, accuracy))
	print("Estimated memory requirement: {0:.2f} MBytes".format(2*8*shape.xpoints*shape.ypoints*N/1024.0/1024.0))

	global f_round
	f_circ=np.zeros((shape.xpoints,shape.ypoints),dtype=np.complex128)
	f_round=np.zeros((shape.xpoints,shape.ypoints,N),dtype=np.complex128)
      
	# move impinging field into cavity
	f_circ = np.sqrt(kat.itmX.T.value) * laser
	# this counts as the first (zeroth) roundtrip
	f_round[:,:,0] = f_circ

	print(" --- pre computing all rountrip fields ---")
	# This will take some time, let's show a progress bar
	p = ProgressBar(maxval=N, widgets=["f_circ:", Percentage(),"|", Timer(), "|", ETA(), Bar()])

	for n in range(1,N):
		#f_circ = FFT_propagate(f_circ,shape,Lambda,LX,1) 
		f_circ = FFT_propagate_simple(f_circ,shape.xpoints,shape.ypoints, shape.xstep, shape.ystep,Lambda,LX,1) 
		f_circ = np.sqrt(kat.etmX.R.value)*FFT_apply_map(f_circ, etm, Lambda)
		#f_circ = FFT_propagate(f_circ,shape,Lambda,LX,1) 
		f_circ = FFT_propagate_simple(f_circ,shape.xpoints,shape.ypoints, shape.xstep, shape.ystep,Lambda,LX,1) 
		f_circ = np.sqrt(kat.itmX.R.value)*FFT_apply_map(f_circ, itm, Lambda)
		f_round[:,:,n] = f_circ;
		p.update(n)
	p.finish()

	print(" --- saving data to file ---")
	import time
	timestr = time.strftime("%Y:%m:%d-%H:%M:%S")
	np.save('fround-'+timestr,f_round)
	print("Saved data into file: {0}".format('fround-'+timestr))
	
	
def FFT_apply_map(field, Map, Lambda):
	k=2.0*np.pi/Lambda
	return field*np.exp(-1j * 2.0 * k * Map.data*Map.scaling);

def plot_field(field):
	ax,fig=plot_setup()
	im = ax.imshow(np.abs(field),origin='lower', aspect='auto')
	cb = fig.colorbar(im, format="%.4g")
	#cb.set_clim(-1.0*zlimit, zlimit)
	ax.autoscale(False)#
	plt.show(block=0)

def plot_propagation(field, grid, Lambda, axis, normalise, Lrange, nr):
	# axis = 0 for x and =1 for y crosssection
	# if normalise = 1, each slice is normalises to peak intensity
	[n,m]=field.shape
	slices = np.zeros((n,np.size(Lrange)))
	from pykat.optics.fft import FFT_propagate
	for idx, L in enumerate(Lrange):
		field2 = FFT_propagate(field,grid, Lambda, L, nr)
		if axis==0:
			slices[:,idx]=np.abs(field2[:,int(m/2)])**2
		else:
			slices[:,idx]=np.abs(field2[int(n/2),:])**2
		if normalise==1:
			peak_intensity= np.abs(field2[int(n/2),int(m/2)])**2
			slices[:,idx]=slices[:,idx]/peak_intensity
	ax,fig=plot_setup()
	im = ax.imshow(slices, aspect='auto')
	cb = fig.colorbar(im, format="%.4g")
	#cb.set_clim(-1.0*zlimit, zlimit)
	ax.autoscale(False)
	plt.show(block=0)


if __name__ == '__main__':
    main()

