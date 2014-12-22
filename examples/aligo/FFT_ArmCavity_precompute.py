from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

import copy
from collections import namedtuple
from collections import OrderedDict
import shelve

import pylab as pl
from  pykat.utilities.plotting.tools import printPDF
from pykat.external.progressbar import ProgressBar, ETA, Percentage, Bar


from pykat.optics.maps import *
from pykat.optics.gaussian_beams import HG_beam, beam_param
from pykat.optics.fft import *

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
	Lambda = aligo.Lambda

	# load and create mirror maps
	global itm, etm
	surface=read_map('etm08_virtual.txt')
	itm=curvedmap('itm_Rc',surface.size,surface.step_size, aligo.itmX_Rc)
	etm=curvedmap('etm_Rc',surface.size,surface.step_size, aligo.etmX_Rc)

	# apply measured map to etm
	#etm.data = etm.data + surface.data

	# setup grid for FFT propagation
	[xpoints,ypoints] = surface.size
	xsize = xpoints * surface.step_size[0]
	ysize = 0.05 * surface.step_size[0]
	xoffset = 0.0
	yoffset = 0.0

	global shape
	shape = grid(xpoints, ypoints, xsize, ysize, xoffset, yoffset)
	x = shape.xaxis
	y = shape.yaxis

	result['shape']=shape
	
	# generate roughly mode-matched input beam
	gx = beam_param(w0=0.012, z=-1834.0)
	gy = gx
	beam = HG_beam(gx,gy,0,0)
	laser = beam.Unm(x,y) 

	R=aligo.etmX_R*aligo.itmX_R
	Loss = 1-R
	accuracy=100E-6
	print("cavity loss: {0}".format(Loss))	
	N=int(required_roundtrips(Loss,accuracy))
	print("required rountrips: {0} (for accuracy of {1})".format(N, accuracy))
	print("Estimated memory requirement: {0:.2f} MBytes".format(2*8*xpoints*ypoints*N/1024.0/1024.0))

	global f_round
	f_circ=np.zeros((xpoints,ypoints),dtype=np.complex128)
	f_round=np.zeros((xpoints,ypoints,N),dtype=np.complex128)
      
	# move impinging field into cavity
	f_circ = np.sqrt(aligo.itmX_T) * laser
	# this counts as the first (zeroth) roundtrip
	f_round[:,:,1] = f_circ

	print(" --- pre computing all rountrip fields ---")
	# This will take some time, let's show a progress bar
	p = ProgressBar(maxval=N, widgets=["computing f_circ:", Percentage(),"|", ETA(), Bar()])

	for n in range(2,N):
		f_circ = FFT_propagate(f_circ,shape,Lambda,aligo.LX,1) 
		f_circ = aligo.etmX_r*FFT_apply_map(f_circ, etm, Lambda)
		f_circ = FFT_propagate(f_circ,shape,Lambda,aligo.LX,1) 
		f_circ = aligo.itmX_r*FFT_apply_map(f_circ, itm, Lambda)
		f_round[:,:,n] = f_circ;
		p.update(n)

	print(" --- saving data to file ---")
	import time
	timestr = time.strftime("%Y:%m:%d-%H:%M:%S")
	np.save('fround-'+timestr,f_round)

	# now the result variables:
	tmpfile = shelve.open(tmpresultfile)
	tmpfile['result']=result
	tmpfile.close()
	
	
def FFT_apply_map(field, Map, Lambda):
	k=2.0*np.pi/Lambda
	return field*np.exp(-1j * k * Map.data*Map.scaling);

if __name__ == '__main__':
    main()

