"""
------------------------------------------------------
Functions related to FFT propogation of beams.
Work in progress, currently these functions are
untested!

Andreas 30.11.2014
http://www.gwoptics.org/pykat/
------------------------------------------------------
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np
import math

def apply_lens(field, grid, Lambda, f):
	# apply a phase factor representing a lens
	k= 2.0*np.pi/Lambda
	return field*(np.exp(2.0 * 1j * k * (2*f - np.sign(f)*np.sqrt((2.0*f)**2-grid.r_squared))))

def apply_thin_lens(field, grid, Lambda, f):
	# apply a phase factor representing a thin lens
	k= 2.0*np.pi/Lambda
	return field*(np.exp(1.0 * 1j * k * grid.r_squared/(2.0*f)))
	
def FFT_propagate(field, grid, Lambda, distance, nr):
	# FFT propagation code in a fixed grid

	k = 2.0*np.pi/Lambda*nr
	plD = np.pi*Lambda*distance/nr

	field = np.fft.fft2(field)
	field = field * np.exp(-1j*k*distance) * np.exp(1j*plD*grid.fft_ir_squared)
	field = np.fft.ifft2(field)

	return field

def FFT_scale_propagate(field, grid0, grid1, Lambda, distance, w0, w1, nr):
	# FFT propagation code with an adaptive grid size.
	# Propagates to a scaled coordinate system, see Virgo Book of
	# Physics pages 179-184, the scaling factor is given
	# as w1/w0 with w0 the beam size at the start of propagation
	# and w1 the expected beam size at the end of propatation.
	# NOT YET TESTED

	k = 2.0*np.pi/Lambda*nr
	plD = np.pi*Lambda*distance*w0/w1/nr
	z0 = distance/(w1/w0-1.0)
	
	# initial scaling
	field = field * np.exp(-1j*k*grid0.r_squared/(2.0*z0))
	field = np.fft.fft2(field)
	# scaled propagator
	field = field * np.exp(-1j*k*distance) * np.exp(1j*plD*grid0.fft_ir_squared)
	field = np.fft.ifft2(field)
	# final scaling
	field = field *w0/w1 * np.exp(1j* grid1.r_squared*(z0+L)/(2.0*z0*z0))
	
	return field

def required_roundtrips(Loss, accuracy):
	# Estimates the required number of roundtrips
	# for a cavity FFT simulation:
	# After Inf roundtrips the power is given by
	# P_total=P_in/(1-r)^2
	# after N roundtrips it is
	# P_N=P_in (1-r^{N+1))^2/(1-r)^2
	# with P_in the power transmitted into the cavity
	# So the relative error is
	# 1-P_N/P_total=1-(1-r^(N+1))^2 which should be smaller than accuracy.
	# This yields
	# roundtrips = ceil(log(1-sqrt(1-accuracy))/log(sqrt(R)))-1;
	# for accuracy<<1 we can simplify this and set 
	R=1-Loss
	return 2*math.ceil(np.log(0.5*accuracy)/np.log(R))
   
class grid():
	# Data structure to describe the size and axes for a (x,y) data array
	# of complex beam amplitudes. Also contain also data structures for
	# FFT propagation
	
	def __init__ (self,xpoints, ypoints, xsize, ysize, xoffset, yoffset):

		self.xpoints=xpoints
		self.ypoints=ypoints
		self.xsize=xsize
		self.ysize=ysize
		self.xoffset=xoffset
		self.yoffset=yoffset

		# compute x and y axis
		self.xstep=self.xsize/self.xpoints
		self.ystep=self.ysize/self.ypoints
		xvector= np.arange(self.xpoints)
		yvector= np.arange(self.ypoints)
		self.xaxis=-self.xsize/2.0 + self.xstep/2.0 + xvector*self.xstep + self.xoffset
		self.yaxis=-self.ysize/2.0 + self.ystep/2.0 + yvector*self.ystep + self.yoffset

		# and some useful variables based on the axis
		self.X,self.Y = np.meshgrid(self.xaxis,self.yaxis)
		self.r_squared = (self.X)**2 + (self.Y)**2
		self.r = np.sqrt(self.r_squared)
		self.angle = np.arctan2(self.Y,self.X)

		# compute frequency axis
		self.xaxis_fft = np.fft.fftshift(np.fft.fftfreq(self.xpoints))/self.xstep
		self.yaxis_fft = np.fft.fftshift(np.fft.fftfreq(self.ypoints))/self.ystep

		# some useful variables based on the frequency axis
		self.fft_X,self.fft_Y = np.meshgrid(self.xaxis_fft, self.yaxis_fft)
		self.fft_ir_squared= np.fft.ifftshift((self.fft_X)**2+(self.fft_Y)**2)
