from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from pykat import finesse
from pykat.commands import *
import pylab as pl
import scipy
from scipy.optimize import fmin
import numpy as np
import pickle
import copy
import sys


def main():
	print("""
	--------------------------------------------------------------
	Example file for using PyKat to automate Finesse simulations
	Finesse: http://www.gwoptics.org/finesse
	PyKat:	 http://www.gwoptics.org/pykat

	The file runs through the various Finesse simulations
	to generate the Finesse results reported in the document:
	`Comparing Finesse simulations, analytical solutions and OSCAR 
	simulations of Fabry-Perot alignment signals', LIGO-T1300345,
	freely available online: http://arxiv.org/abs/1401.5727

	This file is part of a collection; it outputs the results
	shown the document's sections 5 and 6 and saves temporary
	data and a new Finesse input file to be read by master3.py,
	and master4.py.
	
	Andreas Freise 16.01.2014
	--------------------------------------------------------------
	""")
	
	# shall we clear the workspace?
	# %reset -f

	# making these global during testing and debugging
	#global kat, out
	
	kat = finesse.kat(tempdir=".",tempname="test")
	kat.verbose = False
	
	tmpresultfile = "myshelf1.dat"
	
	# loading data saved by master.py
	kat.loadKatFile('asc_base2.kat')
	try:
		with open(tmpresultfile, 'rb') as handle:
			result = pickle.load(handle)
	except: raise Exception("Could not open temprary results file {0}".format(tmpresultfile))
		
	# overwriting some variables
	kat.maxtem=3
	Lambda=1064.0e-9

	# disable PDH photo diode as we won't need it for most of this
	kat.PDrefl_p.enabled = False
	kat.PDrefl_q.enabled = False

	# simulating a tuned cavity
	kat.ETM.phi=result['phi_tuned']
	
	print("--------------------------------------------------------")
	print(" 5. checking wavefronts for ITM/ETM tilt of 0.1nrad")
	tilt(kat)
	
	print("--------------------------------------------------------")
	print(" 6. compute beam tilt from beam propogation")
	gravity_tilt(kat)

	print("--------------------------------------------------------")
	print(" 7. compute optimal demodulation phase of WFS1 and WFS2")

	# adding wave front sensors to global kat object, will need them later
	# on as well.
	
	code_WFS1 = """
	pd1 WFS1_I 9M 0 nWFS1
	pdtype WFS1_I y-split
	pd1 WFS1_Q 9M 90 nWFS1
	pdtype WFS1_Q y-split
	scale 2 WFS1_I % compensate the 0.5 gain of the demodulation
	scale 2 WFS1_Q % compensate the 0.5 gain of the demodulation
	"""
	code_WFS2 = """
	pd1 WFS2_I 9M 0 nWFS2
	pdtype WFS2_I y-split
	pd1 WFS2_Q 9M 90 nWFS2
	pdtype WFS2_Q y-split
	scale 2 WFS2_I % compensate the 0.5 gain of the demodulation
	scale 2 WFS2_Q % compensate the 0.5 gain of the demodulation
	"""
	kat.parseKatCode(code_WFS1)
	kat.parseKatCode(code_WFS2)
	
	(WFS1_phase, WFS2_phase) = asc_phases(kat)
	kat.WFS1_I.phase1=WFS1_phase
	kat.WFS1_Q.phase1=WFS1_phase+90.0
	kat.WFS2_I.phase1=WFS2_phase
	kat.WFS2_Q.phase1=WFS2_phase+90.0
	result['WFS1_phase']=WFS1_phase
	result['WFS2_phase']=WFS2_phase

	print("--------------------------------------------------------")
	print(" 8. compute ASC signal matrix at WFS1 and WFS2")
	signal = asc_signal(kat)
	
	print("--------------------------------------------------------")
	print(" Saving results in temp. files to be read by master3.py")
	# re-enable PDH photo diode for saving of kat file for next script
	kat.PDrefl_p.enabled = True
	kat.PDrefl_q.enabled = True

	tmpkatfile = "asc_base3.kat"
	tmpresultfile = "myshelf2.dat"
	print(" kat object saved in: {0}".format(tmpkatfile))
	print(" current results saved in: {0}".format(tmpresultfile))
	kat.saveScript(tmpkatfile)
	with open(tmpresultfile, 'wb') as handle:
		pickle.dump(result, handle)


#-----------------------------------------------------------------------------------

def asc_signal(tmpkat):
	kat = copy.deepcopy(tmpkat)

	code_lock = """
	set err PDrefl_p re
	lock z $err 900 1p
	put* ETM phi $z
	noplot z
	"""
	
	kat.parseKatCode(code_lock)
	# need to re-enable the photo diode for lock
	kat.PDrefl_p.enabled = True

	kat.parseKatCode('yaxis abs')
	kat.noxaxis = True
	kat.maxtem=1

	signal=np.zeros((2, 2))
	kat.ITM.ybeta=1e-10
	kat.ETM.ybeta=0.0
	out = kat.run()
	signal[0,0] = out["WFS1_I"]
	signal[1,0] = out["WFS2_I"]

	
	kat.ITM.ybeta=0.0
	kat.ETM.ybeta=-1e-10
	out = kat.run()
	signal[0,1] = out["WFS1_I"]
	signal[1,1] = out["WFS2_I"]
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
		kat.WFS1_I.phase1=x[0]
		out = kat.run()
		signal = out["WFS1_I"]
		print('\r minimising: function value {0:<16g}'.format(float(signal)), end='')
		sys.stdout.flush()
		return -1*abs(signal)

	def demod_phase2(x):
		kat.WFS2_I.phase1=x[0]
		out = kat.run()
		signal = out["WFS2_I"]
		print('\r minimising: function value {0:<16g}'.format(float(signal)), end='')
		sys.stdout.flush()
		return -1*abs(signal)

	kat.ITM.ybeta=1e-10
	kat.ETM.ybeta=0.0
	res = fmin(demod_phase1, [0.0], xtol=1e-8, disp=False)
	WFS1_phase = res[0]
	print("")
	print(" WFS1 demod phase : %.10g deg" % WFS1_phase)
	 
	kat.ITM.ybeta=0.0
	kat.ETM.ybeta=-1e-10
	res = fmin(demod_phase2, [0.0], xtol=1e-8, disp=False)
	WFS2_phase = res[0]
	print("")
	print(" WFS2 demod phase : %.10g deg" % WFS2_phase)
	return(WFS1_phase, WFS2_phase)	  
	
def gravity_tilt(tmpkat):
	kat = copy.deepcopy(tmpkat)

	def compute_gravity_tilt(tmpkat):
		#global out, y1, y2
		kat = copy.deepcopy(tmpkat)
		out = kat.run()
		
		y1 = np.abs(out["b1"])
		y2 = np.abs(out["b1_1k"])
		# shift of beam center	on detector 1 (as m/w0y)
		x1 = np.sum(out.x*y1)/np.sum(y1) 
		# shift of beam center	on detector 2 (as m/w0y)
		x2 = np.sum(out.x*y2)/np.sum(y2)
		# calibrate this in meter by mutliplying with w0y
		# and compute the angle geometrically		 
		w0=out["w0y"][0]
		detector_distance = 1000.0
		tilt=w0*(x2-x1)/detector_distance
		print(" Wavefront tilt : %g rad" % tilt)

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
	print(" ITM ybeta 0.1nrad")
	kat.parseKatCode(code_WFS1)
	kat.parseKatCode(code_xaxis)
	kat.spo1.L=1000.0
	kat.ITM.ybeta=1e-10
	kat.ETM.ybeta=0.0
	compute_gravity_tilt(kat)
	print(" ETM ybeta -0.1nrad")
	kat.ITM.ybeta=0.0
	kat.ETM.ybeta=-1e-10
	compute_gravity_tilt(kat)

	print(" WFS2:")
	print(" ITM ybeta 0.1nrad")
	kat = copy.deepcopy(tmpkat)
	kat.parseKatCode(code_WFS2)
	kat.parseKatCode(code_xaxis)
	kat.spo1.L=1.0e-9
	kat.ITM.ybeta=1e-10
	kat.ETM.ybeta=0.0
	compute_gravity_tilt(kat)
	print(" ETM ybeta -0.1nrad")
	kat.ITM.ybeta=0.0
	kat.ETM.ybeta=-1e-10
	compute_gravity_tilt(kat)

	
def tilt(tmpkat):
	kat = copy.deepcopy(tmpkat)
	
	def compute_tilt(tmpkat):
		#global kat, out
		kat = copy.deepcopy(tmpkat)
		out = kat.run()
		# compute data x range in meters
		beamsize = out["w0y"][0].real
		xrange = beamsize*(out.x.max()-out.x.min())
		stepsize=xrange/(len(out.x)-1)
		print(" Beamsize %e m" % beamsize)
		print(" Measurement range: %e m, stepszie: %e m" % (xrange, stepsize))
		# compute difference in angle between wavefront of carrier and sidebands
		#global diff_l, diff_u, tilt_l, tilt_u
		diff_l = (np.angle(out["PDrefl_low"], deg=True) - np.angle(out["PDrefl_car"], deg=True))/stepsize
		diff_u = (np.angle(out["PDrefl_up"], deg=True)  - np.angle(out["PDrefl_car"], deg=True))/stepsize
		tilt_l = diff_l[1:-1]-diff_l[0:-2]
		tilt_u = diff_u[1:-1]-diff_u[0:-2]
		print(" Tilt (upper	 - car), mean: %e m/deg, stddev %e m/deg" % (np.mean(tilt_u), np.std(tilt_u)))
		print(" Tilt (lower	 - car), mean: %e m/deg, stddev %e m/deg" % (np.mean(tilt_l), np.std(tilt_l)))
		return (np.mean(tilt_l), np.mean(tilt_u))

	code_WFS1 = """
	beam PDrefl_car 0 nWFS1
	beam PDrefl_up 9M nWFS1
	beam PDrefl_low -9M nWFS1
	bp w0y y w0 nWFS1
	"""

	code_WFS2 = """
	beam PDrefl_car 0 nWFS2
	beam PDrefl_up 9M nWFS2
	beam PDrefl_low -9M nWFS2
	bp w0y y w0 nWFS2
	"""
	code_comm = """
	xaxis PDrefl_car y lin -1 1 100
	put PDrefl_up y $x1
	put PDrefl_low y $x1
	yaxis abs:deg
	"""

	print(" WFS1:")
	print(" ITM ybeta 0.1nrad")
	kat.parseKatCode(code_comm)
	kat.parseKatCode(code_WFS1)
	kat.ITM.ybeta=1e-10
	kat.ETM.ybeta=0.0
	(a1, a2) = compute_tilt(kat)
	
	print(" ETM ybeta -0.1nrad")
	kat.ITM.ybeta=0.0
	kat.ETM.ybeta=-1e-10
	(a3, a4) = compute_tilt(kat)
	
	print(" WFS2:")
	print(" ITM ybeta 0.1nrad")
	kat = copy.deepcopy(tmpkat)
	kat.parseKatCode(code_comm)
	kat.parseKatCode(code_WFS2)
	kat.ITM.ybeta=1e-10
	kat.ETM.ybeta=0.0
	(a5, a6) = compute_tilt(kat)

	print(" ETM ybeta -0.1nrad")
	kat.ITM.ybeta=0.0
	kat.ETM.ybeta=-1e-10
	(a6, a7) = compute_tilt(kat)

	return 
	
if __name__ == '__main__':
	main()

