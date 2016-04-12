from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from pykat import finesse
from pykat.commands import *
from pykat.optics.gaussian_beams import beam_param

import pylab as pl
import scipy
from scipy.optimize import minimize_scalar
import numpy as np
import shelve
import copy
import sys

# making python2 and python3 compatible
try:
   input = raw_input
except NameError:
   pass

def main():

	print("""
	--------------------------------------------------------------
	Example file for using PyKat to automate Finesse simulations
	Finesse: http://www.gwoptics.org/finesse
	PyKat:	 https://pypi.python.org/pypi/PyKat/
	
	The file runs through the various pykat files which are used
	to generate the Finesse results reported in the document:
	`Comparing Finesse simulations, analytical solutions and OSCAR 
	simulations of Fabry-Perot alignment signals', LIGO-T1300345
	
	This file is part of a collection. Run this after master2.py
	
	Andreas Freise 06.12.2013
	--------------------------------------------------------------
	""")
	
	# shall we clear the workspace?
	# %reset -f

	# making these global during testing and debugging
	global kat
	global out

	kat = finesse.kat(tempdir=".",tempname="test")
	kat.verbose = False

	tmpresultfile = 'myshelf2.dat'
	
	# loading data saved by master.py
	kat.loadKatFile('asc_base3.kat')
	try:
		tmpfile = shelve.open(tmpresultfile)
		result=tmpfile[str('result')]
		tmpfile.close()
	except: raise Exception("Could not open temprary results file {0}".format(tmpresultfile))

		
	# overwriting some variables
	kat.maxtem=3
	Lambda=1064.0e-9
	
	# this does not work yet due to the scale command
	kat.PDrefl_p.enabled = False
	kat.PDrefl_q.enabled = False
	kat.WFS1_I.enabled = False
	kat.WFS1_Q.enabled = False
	kat.WFS2_I.enabled = False
	kat.WFS2_Q.enabled = False

	kat.ETM.phi=result['phi_tuned']

	(beam1, beam2, beam3) = get_qs(kat)
	"""
	print "	 Measured beam parameter:" 
	print "	 - At front of ITM (no thermal lens) q={0}".format(beam1.q)
	print "	   (eqals w0={0}, z={1})".format(beam1.w0, beam1.z)
	print "	 - At pick of mirror 'po' (50k lens) q={0}".format(beam2.q)
	print "	   (eqals w0={0}, z={1})".format(beam2.w0, beam2.z)
	print "	 - At pick of mirror 'po' (5k lens) q={0}".format(beam3.q)
	print "	   (eqals w0={0}, z={1})".format(beam3.w0, beam3.z)
	#print "  Setting these now view Gauss command and adding thermal lens"
	"""
	kat.ITM.nITM1.node.setGauss(kat.ITM,beam1)

	print("--------------------------------------------------------")
	print(" 11. computing beam sizes  with thermal lens")
	#beam_size(kat, beam2, beam3)
	
	kat.ITM_TL.f=50e3
	kat.maxtem = 8
	print("--------------------------------------------------------")
	print(" 11. computing beam tilt with thermal lens (f={0}, maxtem={1})".format(kat.ITM_TL.f, kat.maxtem))
	#gravity_tilt(kat)

	kat.ITM_TL.f=5e3
	kat.maxtem = 23

	print("--------------------------------------------------------")
	print(" 12. computing beam tilt with thermal lens (f={0}, maxtem={1})".format(kat.ITM_TL.f, kat.maxtem))
	#gravity_tilt(kat)

	
	print("--------------------------------------------------------")
	print(" 12. compute beam center with thermal lens")

	

	
	print("--------------------------------------------------------")
	print(" Saving results in temp. files to be read by master6.py")
	tmpkatfile = "asc_base4.kat"
	tmpresultfile = "myshelf3.dat"
	print(" kat object saved in: {0}".format(tmpkatfile))
	print(" current results saved in: {0}".format(tmpresultfile))
	# first the current kat file
	kat.saveScript(tmpkatfile)
	# now the result variables:
	tmpfile = shelve.open(tmpresultfile)
	tmpfile[str('result')]=result
	tmpfile.close()


#-----------------------------------------------------------------------------------

def get_qs(tmpkat):
	kat = copy.deepcopy(tmpkat)
	nodename0="npsl"
	nodename1="nITM1"
	nodename2="nWFS1"
	nodename3="nWFS2"
	# measure beam parameter for the 'cold beam' i.e. the laser beam
	# matched to the cavity without any thermal lens
	code_bp = "bp w0 y q {0}\nbp w1 y q {1}\nbp w2 y q {2}\nbp w3 y q {3}".format(nodename0,nodename1,nodename2,nodename3)

	kat.parseKatCode(code_bp)
	kat.parseKatCode('yaxis re:im')
	kat.noxaxis = True
	kat.maxtem=0

	def beam_size(tmpkat, f):
		kat = copy.deepcopy(tmpkat)
		# 1. run finesse with input laser mode matched to cavity (no thermal lens)
		out = kat.run()

		# beam at laser when matched to cold cavity
		# (note the sign flip of the real part to change direction of gauss param)
		q0 = -1.0*out['w0'].conjugate()
		beam0 = beam_param(q=q0)
		kat.psl.npsl.node.setGauss(kat.psl, beam0)
		kat.parseKatCode("startnode npsl")

		# add thermal lens and propagate input beam to ITM
		kat.ITM_TL.f=f
		if "ITM_TL_r" in kat._kat__components:
			kat.ITM_TL_r.f=f
		out = kat.run()
		
		# computing beam size at ITM 
		# and then we reflect of ITM, an set it as new startnode
		q_in = out['w1']
		from pykat.optics.ABCD import apply, mirror_refl
		abcd = mirror_refl(1,-2500)
		q_out = apply(abcd,q_in,1,1)
		beam1 = beam_param(q=q_out)	
		kat.removeLine("startnode")
		kat.psl.npsl.node.removeGauss()
		if "ITM_TL_r" in kat._kat__components:
			kat.ITM.nITM1r.node.setGauss(kat.ITM, beam1)
			kat.parseKatCode("startnode nITM1r")
		else:
			kat.ITM.nITM1.node.setGauss(kat.ITM, beam1)
			kat.parseKatCode("startnode nITM1")
		out = kat.run()

		# computing beam size at WFS1 and WFS2
		q2 = out['w2']
		beam2 = beam_param(q=q2)	 
		q3 = out['w3']
		beam3 = beam_param(q=q3)	 
		print("	 Sideband (input mode) beam size with thermal lens f={0}".format(f))
		print("	 - WFS1 w={0:.6}cm".format(100.0*beam2.w))
		print("	   (w0={0}, z={1})".format(beam2.w0, beam2.z))
		print("	 - WFS2 w={0:.6}cm".format(100.0*beam3.w))
		print("	   (w0={0}, z={1})".format(beam3.w0, beam3.z))
		input("Press enter to continue")
		return(beam1, beam2, beam3)
		
	f=50e3
	beam_size(kat,f)
	f=5e3
	(beam1,beam2,beam3)=beam_size(kat,f)
	return (beam1, beam2,beam3)

def asc_signal(tmpkat):
	kat = copy.deepcopy(tmpkat)

	code_lock = """
	set err PDrefl_p re
	lock z $err 900 1p
	put* ETM phi $z
	noplot z
	"""
	
	kat.parseKatCode(code_lock)
	kat.parseKatCode('yaxis abs')
	kat.noxaxis = True
	kat.maxtem=1

	signal=np.zeros((2, 2))
	kat.ITM.ybeta=1e-10
	kat.ETM.ybeta=0.0
	out = kat.run()
	WFS1_idx=out.ylabels.index("WFS1_I")
	WFS2_idx=out.ylabels.index("WFS2_I")
	signal[0,0] = out.y[WFS1_idx]
	signal[1,0] = out.y[WFS2_idx]

	kat.ITM.ybeta=0.0
	kat.ETM.ybeta=-1e-10
	out = kat.run()
	signal[0,1] = out.y[WFS1_idx]
	signal[1,1] = out.y[WFS2_idx]
	signal = signal *1e10
	sensors=('WFS1', 'WFS2')
	mirrors=('ITM', 'ETM')
	print("	 ASC Matrix:")
	for i in range(2):
		print("	 ", sensors[i], " ", end=' ')
		for j in range(2):
			print("%12.10g" % signal[i,j], end=' ')
		print(mirrors[i])
	return signal
	
	
def asc_phases(tmpkat):
	kat = copy.deepcopy(tmpkat)
	
	kat.parseKatCode('yaxis abs')
	kat.noxaxis = True
	kat.maxtem=1

	def demod_phase1(x):
		kat.WFS1_I.phi[0]=x
		out = kat.run()
		WFS1_idx=out.ylabels.index("WFS1_I")
		signal = out.y[WFS1_idx]
		print('\r minimising: function value %g					   ' % signal, end=' ')
		sys.stdout.flush()
		return -1*abs(signal)

	def demod_phase2(x):
		kat.WFS2_I.phi[0]=x
		out = kat.run()
		WFS2_idx=out.ylabels.index("WFS2_I")
		signal = out.y[WFS2_idx]
		print('\r minimising: function value %g					   ' % signal, end=' ')
		sys.stdout.flush()
		return -1*abs(signal)

	kat.ITM.ybeta=1e-10
	kat.ETM.ybeta=0.0
	res = minimize_scalar(demod_phase1, method='brent')
	WFS1_phase = res.x
	print("")
	print(" WFS1 demod phase : %.10g deg" % WFS1_phase)
	 
	kat.ITM.ybeta=0.0
	kat.ETM.ybeta=-1e-10
	res = minimize_scalar(demod_phase2, method='brent')
	WFS2_phase = res.x
	print("")
	print(" WFS2 demod phase : %.10g deg" % WFS2_phase)
	return(WFS1_phase, WFS2_phase)	  
	
def gravity_tilt(tmpkat):
	kat = copy.deepcopy(tmpkat)

	def compute_gravity_tilt(tmpkat):
		kat = copy.deepcopy(tmpkat)
		out = kat.run()

		y1 = out["b1"]
		y2 = out["b1_1k"]
		# shift of beam center	on detector 1 (as m/w0y)
		x1 = np.sum(out.x*y1)/np.sum(y1) 
		# shift of beam center	on detector 2 (as m/w0y)
		x2 = np.sum(out.x*y2)/np.sum(y2)
		# calibrate this in meter by mutliplying with w0y
		# and compute the angle geometrically		 
		w0=out["w0y"][0]
		detector_distance = 1000.0
		tilt=w0*(x2-x1)/detector_distance
		print(" Wavefront tilt : %g nrad" % tilt)

	code_WFS1 = """
	beam b1 nWFS1
	beam b1_1k nL1_in
	bp w0y y w0 nWFS1
	"""

	code_WFS2 = """
	m md 0 1 0 nWFS2 nWFS2b
	s sd 1k nWFS2b nWFS2c
	beam b1 nWFS2*
	beam b1_1k nWFS2c
	bp w0y y w0 nWFS2
	"""

	code_xaxis= """
	xaxis b1 y lin -40 40 800
	put b1_1k y $x1
	yaxis abs
	"""
	print(" WFS1:")
	print(" ITM ybeta 0.1nm")
	kat.parseKatCode(code_WFS1)
	kat.parseKatCode(code_xaxis)
	kat.spo1.L=1000.0
	kat.ITM.ybeta=1e-10
	kat.ETM.ybeta=0.0
	compute_gravity_tilt(kat)
	print(" ETM ybeta -0.1nm")
	kat.ITM.ybeta=0.0
	kat.ETM.ybeta=-1e-10
	compute_gravity_tilt(kat)

	print(" WFS1:")
	print(" ITM ybeta 0.1nm")
	kat = copy.deepcopy(tmpkat)
	kat.parseKatCode(code_WFS2)
	kat.parseKatCode(code_xaxis)
	kat.spo1.L=1.0e-9
	kat.ITM.ybeta=1e-10
	kat.ETM.ybeta=0.0
	compute_gravity_tilt(kat)
	print(" ETM ybeta -0.1nm")
	kat.ITM.ybeta=0.0
	kat.ETM.ybeta=-1e-10
	compute_gravity_tilt(kat)

def beam_size(tmpkat, beam2, beam3):
	kat = copy.deepcopy(tmpkat)

	global out
	code_bps = """
	bp wWFS1 y w nWFS1
	bp wWFS2 y w nWFS2
	"""
	kat.parseKatCode(code_bps)
	kat.maxtem = 0
	kat.ITM.R=1.0 
	kat.ITM.T=0.0 
	
	kat.noxaxis = True
	kat.ITM_TL.f=50e3
	if "ITM_TL_r" in kat._kat__components:
		kat.ITM_TL_r.f=50e3

	kat.po.nWFS1.node.setGauss(kat.po,beam2)

	out = kat.run()

	WFS1_idx=out.ylabels.index("wWFS1")
	WFS2_idx=out.ylabels.index("wWFS2")
	
	
	y1 = out.y[WFS1_idx]
	y2 = out.y[WFS2_idx]
	print("	 Beam size with thermal lens f={0}".format(kat.ITM_TL.f))
	print("	 WFS1: {0}cm".format(y1*100.0))
	print("	 WFS2: {0}cm".format(y2*100.0))
	
	kat.ITM_TL.f=5e3
	if "ITM_TL_r" in kat._kat__components:
		kat.ITM_TL_r.f=5e3
	kat.po.nWFS1.node.setGauss(kat.po,beam3)
	out = kat.run()
	y1 = out.y[WFS1_idx]
	y2 = out.y[WFS2_idx]
	print("	 Beam size with thermal lens f={0}".format(kat.ITM_TL.f))
	print("	 WFS1: {0}cm".format(y1*100.0))
	print("	 WFS2: {0}cm".format(y2*100.0))

	
if __name__ == '__main__':
	main()

