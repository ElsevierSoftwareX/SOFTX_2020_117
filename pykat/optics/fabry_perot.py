# ------------------------------------------------------
# Moving some utility function from Simtools
# http://www.gwoptics.org/simtools/
# to PyKat. Work in progress, currently these
# functions are not tested.
# 
# Andreas 30.11.2014
# ------------------------------------------------------
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


import numpy as np

def cavity_info(Lambda, L, Rc1, Rc2):
	g1 = 1-L/Rc1
	g2 = 1-L/Rc2
	G  = g1*g2
	G1 = (g1+g2-2*G)
	
	k = 2.0*np.pi/Lambda

	zr = np.sqrt( L**2 * G * (1.0-G)/G1**2) 
	z1 = L * g2 * (1-g1) / G1
	z2 = L * g1 * (1-g2) / G1
	w0 = L * Lambda / np.pi * np.sqrt( G * (1-G) / G1**2 )
	w1 = np.sqrt( L * Lambda /np.pi * np.sqrt( g2 / (g1 * (1-G)) ))
	w2 = np.sqrt( L * Lambda /np.pi * np.sqrt( g1 / (g2 * (1-G)) ))
	return [zr, w0, z1, w1, w2]
	
def cavity_w1w2_Rc1Rc2(Lambda, L, w1, w2):
	C = (Lambda*L/np.pi)**2
	g1 = -1.0 * np.sqrt(1-C/(w1*w2)**2)*(w2/w1)
	g2 = -1.0 * np.sqrt(1-C/(w1*w2)**2)*(w1/w2)
	Rc1 = L/(1-g1)
	Rc2 = L/(1-g2)
	return [g1, g2, Rc1, Rc2]
