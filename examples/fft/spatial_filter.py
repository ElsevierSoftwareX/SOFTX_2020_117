from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pylab as pl

from pykat.optics.gaussian_beams import HG_beam, beam_param
from pykat.optics.fft import *
import numpy as np
import scipy

def main():
	# wavelength
	Lambda = 1064.0E-9 
	# distance to propagate/focal length of lens
	D = 10
	# mix coefficients
	c1 = 0.7
	c2 = 0.3
	# mode indices 
	n1 = 2
	m1 = 3
	n2 = 1
	m2 = 0

	######## Generate Grid stucture required for FFT propagation ####
	xpoints = 512
	ypoints = 512
	xsize = 0.05
	ysize = 0.05
	# Apply offset such that the center of the beam lies in the
	# center of a grid tile
	xoffset = -0.5*xsize/xpoints
	yoffset = -0.5*ysize/ypoints

	shape = grid(xpoints, ypoints, xsize, ysize, xoffset, yoffset)
	x = shape.xaxis
	y = shape.yaxis

	######## Generates a mixture of fields ################
	gx = beam_param(w0=2e-3, z=0)
	gy = beam_param(w0=2e-3, z=0)
	beam = HG_beam(gx,gy,n1,m1)
	field1 = beam.Unm(x,y) 
	beam2 = HG_beam(gx,gy,n2,m2)
	field2 = beam2.Unm(x,y) 
	global field, laser1, laser2
	field = np.sqrt(c1)*field1 + np.sqrt(c2)*field2
 
	####### Apply phase plate #######################################

	laser1 = field*(np.conjugate(field1))
	laser2 = field*(np.conjugate(field2))

	####### Propagates the field by FFT ##############################
	laser1 = FFT_propagate(laser1,shape,Lambda,D,1) 
	laser2 = FFT_propagate(laser2,shape,Lambda,D,1) 

	f=D
	#laser1 = apply_lens(laser1, shape, Lambda, f) 
	#laser2 = apply_lens(laser2, shape, Lambda, f) 
	laser1 = apply_thin_lens(laser1, shape, Lambda, f) 
	laser2 = apply_thin_lens(laser2, shape, Lambda, f) 

	laser1 = FFT_propagate(laser1,shape,Lambda,D,1) 
	laser2 = FFT_propagate(laser2,shape,Lambda,D,1)
	
	# midpoint computation for even number of points only!
	midx=(xpoints)//2
	midy=(ypoints)//2
	coef1 = np.abs(laser1[midx,midy])  
	coef2 = np.abs(laser2[midx,midy])

	ratio = (coef1/coef2)**2
	pc2 = 1/(1+ratio)
	pc1 = pc2*ratio

	print("c1 {0}, coef1 {1}, error {3} (raw output {2})".format(c1, pc1, coef1, np.abs(c1-pc1)))
	print("c2 {0}, coef2 {1}, error {3} (raw output {2})".format(c2, pc2, coef2, np.abs(c2-pc2)))

	# plot hand tuned for certain ranges and sizes, not automtically scaled
	fig=pl.figure(110)
	fig.clear()
	off1=xpoints/10
	off2=xpoints/6
	pl.subplot(1, 3, 1)
	pl.imshow(abs(field))
	pl.xlim(midx-off1,midx+off1)
	pl.ylim(midy-off1,midy+off1)
	pl.draw()
	pl.subplot(1, 3, 2)
	pl.imshow(abs(laser1))
	pl.xlim(midx-off2,midx+off2)
	pl.ylim(midy-off2,midy+off2)
	pl.draw()
	pl.subplot(1, 3, 3)
	pl.imshow(abs(laser2))
	pl.xlim(midx-off2,midx+off2)
	pl.ylim(midy-off2,midy+off2)
	pl.draw()
	if in_ipython():
		pl.show(block=0)
	else:
		pl.show(block=1)


# testing if the script is run from within ipython
def in_ipython():
    try:
        __IPYTHON__
    except NameError:
        return False
    else:
        return True

if __name__ == '__main__':
	main()

