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
import pickle
import copy
import sys

def set_thermal_lens(kat, f):
	kat.ITM_TL.f=f
	# if a bs-based cavity is used, we need to set second lens
	if "ITM_TL_r" in kat._kat__components:
		kat.ITM_TL_r.f=f
	return (kat)

def main():

	print("""
	--------------------------------------------------------------
	Example file for using PyKat to automate Finesse simulations
	Finesse: http://www.gwoptics.org/finesse
	PyKat:	 https://pypi.python.org/pypi/PyKat/

	The file runs through the various Finesse simulations
	to generate the Finesse results reported in the document:
	`Comparing Finesse simulations, analytical solutions and OSCAR 
	simulations of Fabry-Perot alignment signals', LIGO-T1300345,
	freely available online: http://arxiv.org/abs/1401.5727

	Run this file after master2.py to create data which can be
	plotted using master4_plot.py. Results are saved after 
	each step and plots can be created at any time.
		
	Andreas Freise 16.01.2014
	--------------------------------------------------------------
	""")
	
	# shall we clear the workspace?
	# %reset -f

	# making these global during testing and debugging
	#global kat, out

	kat = finesse.kat(tempdir=".",tempname="test")
	kat.verbose = False

	tmpresultfile = 'myshelf2.dat'
	
	# loading data saved by master.py
	kat.loadKatFile('asc_base3.kat')
	try:
		with open(tmpresultfile, 'rb') as handle:
			result = pickle.load(handle)
	except: raise Exception("Could not open temprary results file {0}".format(tmpresultfile))
	
	# this does not work yet due to the scale command
	kat.PDrefl_p.enabled = False
	kat.PDrefl_q.enabled = False
	kat.WFS1_I.enabled = False
	kat.WFS1_Q.enabled = False
	kat.WFS2_I.enabled = False
	kat.WFS2_Q.enabled = False

	kat.ETM.phi=result['phi_tuned']

	print("--------------------------------------------------------")
	print(" 11. Do beam tracing to measure beam parameters")
	# get beam parameters at nodes: "npsl", "nITM1", "nWFS1", "nWFS2", "npo2"
	global beam1, beam2, beam3
	beam1 = get_qs(kat,1e13)
	beam2 = get_qs(kat,50e3)
	beam3 = get_qs(kat,5e3)

	# starting with cold beam at npsl
	kat.psl.npsl.node.setGauss(kat.psl, beam1[0])
	kat.parseKatCode("startnode npsl")


	# if we use bs-based cavity we try to set good beam
	# parameter for reflected beam, first by
	# computing 'average' beam from cold and hot beams
	x1=0.70
	x2=0.30
	if "ITM_TL_r" in kat._kat__components:
		beam50 = beam_param(z=(x1*beam1[1].z+x2*beam2[1].z), w0=(x1*beam1[1].w0+x2*beam2[1].w0))
		beam5  = beam_param(z=(x1*beam1[1].z+x2*beam3[1].z), w0=(x1*beam1[1].w0+x2*beam3[1].w0))
		node_text = "at ITM->nITM1r"
		t_comp=kat.ITM
		t_node=kat.ITM.nITM1r
	else:
		beam50 = beam_param(z=(x1*beam1[4].z+x2*beam2[4].z), w0=(x1*beam1[4].w0+x2*beam2[4].w0))
		beam5  = beam_param(z=(x1*beam1[4].z+x2*beam3[4].z), w0=(x1*beam1[4].w0+x2*beam3[4].w0))
		node_text = "at s2->npo2"
		t_comp=kat.s2
		t_node=kat.s2.npo2
	
	kat = set_thermal_lens(kat,50.0e3)

	print("--------------------------------------------------------")
	print(" 12. computing beam tilt with thermal lens (f={0})".format(kat.ITM_TL.f))

	print(" Setting compromise beam parameter {0}:\n w0={1}, z={2}".format(node_text, beam50.w0, beam50.z))
	t_node.node.setGauss(t_comp, beam50)
	kat.maxtem=8
	print(" Calculating maxtem = %d " % kat.maxtem)
	tmp = gravity_tilt(kat)

	kat = set_thermal_lens(kat,5e3)
	print("--------------------------------------------------------")
	print(" 13. computing beam tilt with thermal lens (f={0})".format(kat.ITM_TL.f))
	print(" Setting compromise beam parameter {0}:\n w0={1}, z={2}".format(node_text, beam5.w0, beam5.z))
	t_node.node.setGauss(t_comp, beam5)
	#maxtems = [1, 3, 5, 9, 11, 13, 15, 19, 23, 25, 27, 29]
	#maxtems = [1, 3, 5, 9, 11, 13, 15]
	maxtems = [1, 3, 5, 7]
	converge_tilt(kat, maxtems)
   
	kat.PDrefl_p.enabled = True
	kat.WFS1_I.enabled = True
	kat.WFS2_I.enabled = True
	
	kat = set_thermal_lens(kat,50e3)
	print("--------------------------------------------------------")
	print(" 12. computing alignment signal with thermal lens (f={0})".format(kat.ITM_TL.f))
	t_node.node.setGauss(t_comp, beam50)
	#maxtems = [1, 3, 5, 9, 11, 13, 15, 17, 19]
	#maxtems = [1, 3, 5, 9, 11, 13, 15]
	maxtems = [1, 3, 5, 7]
	converge_asc(kat, maxtems, 'asc_signals_50.txt')

	kat = set_thermal_lens(kat,5e3)
	print("--------------------------------------------------------")
	print(" 13. computing alignment signal with thermal lens (f={0})".format(kat.ITM_TL.f))
	t_node.node.setGauss(t_comp, beam5)
	#maxtems = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31]
	#maxtems = [1, 3, 5, 9, 11, 13, 15]
	maxtems = [1, 3, 5, 7]
	converge_asc(kat, maxtems, 'asc_signals_5.txt')

#-----------------------------------------------------------------------------------

def converge_asc(tmpkat, maxtems, filename):
	kat = copy.deepcopy(tmpkat)
	savedata = np.zeros([len(maxtems),5])
	for idx, tem in enumerate(maxtems):
		savedata[idx,0]=tem
		print(" Calculating maxtem = %d " % tem)
		kat.maxtem = tem
		tmp = asc_signal(kat)
		savedata[idx,1:5]=tmp.reshape([1,4])
		print(" Saving results in file: {0}".format(filename))
		np.savetxt(filename, savedata[0:idx+1,:], fmt=b'%.18e', delimiter=' ')	  

def converge_tilt(tmpkat, maxtems):
	kat = copy.deepcopy(tmpkat)
	savedata = np.zeros([len(maxtems),5])
	filename = "thermal_gravity.txt"
	for idx, tem in enumerate(maxtems):
		savedata[idx,0]=tem
		print(" Calculating maxtem = %d " % tem)
		kat.maxtem = tem
		tmp = gravity_tilt(kat)
		savedata[idx,1:5]=tmp
		print(" Saving results in file: {0}".format(filename))
		np.savetxt(filename, savedata[0:idx+1,:], fmt=b'%.18e', delimiter=' ')	  

def get_qs(tmpkat,f):
	global kat, out	
	kat = copy.deepcopy(tmpkat)

	# measure beam parameter for the 'cold beam' i.e. the laser beam
	# matched to the cavity without any thermal lens
	nodenames=["npsl", "nITM1", "nWFS1", "nWFS2", "npo2"]
	for idx, nname in enumerate(nodenames):
		kat.parseKatCode('bp w{0} y q {1}'.format(idx, nname))

	kat.parseKatCode('yaxis re:im')
	kat.noxaxis = True
	kat.maxtem=0
	
	def beam_size(tmpkat, f, beam0):
		kat = copy.deepcopy(tmpkat)
		kat.psl.npsl.node.setGauss(kat.psl, beam0)
		
		# add thermal lens and propagate input beam to ITM
		kat = set_thermal_lens(kat, f)
		out = kat.run()
		# computing beam size at ITM 
		# and then we reflect of ITM, an set it as new startnode
		q_in = out['w1']
		#import pykat.optics.ABCD as ABCD
		#abcd = ABCD.mirror_refl(1,2500)
		#q_out = ABCD.apply(abcd,q_in,1,1)
		beam1 = beam_param(q=q_in)	   
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

		# computing beam size at pick off
		q4 = out['w4']
		beam4 = beam_param(q=q4)	
		print(" Input mode beam size with thermal lens f={0}".format(f))
		print(" - ITM  w={0:.6}cm  (w0={1}, z={2})".format(100.0*beam1.w,beam1.w0, beam1.z))
		print(" - WFS1 w={0:.6}cm  (w0={1}, z={2})".format(100.0*beam2.w,beam2.w0, beam2.z))
		print(" - WFS2 w={0:.6}cm  (w0={1}, z={2})".format(100.0*beam3.w,beam3.w0, beam3.z))
		print(" - npo2 w={0:.6}cm  (w0={1}, z={2})".format(100.0*beam4.w,beam4.w0, beam4.z))
		#raw_input("Press enter to continue")
		return [beam1, beam2, beam3, beam4]

	#global out
	# run finesse with input laser mode matched to cavity (no thermal lens)
	out = kat.run()
	# beam at laser when matched to cold cavity
	# (note the sign flip of the real part to change direction of gauss param)
	q0 = -1.0*out['w0'].conjugate()
	beam0 = beam_param(q=q0)
	# compute beam sizes when tracing this beam back through the system
	(beam1,beam2,beam3, beam4)=beam_size(kat,f,beam0)
	
	return (beam0, beam1, beam2,beam3, beam4)

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
	print(" ASC Matrix:")
	for i in range(2):
		print("	 ", sensors[i], " ", end=' ')
		for j in range(2):
			print("%12.10g" % signal[i,j], end=' ')
		print(mirrors[i])
	return signal
	
def gravity_tilt(tmpkat):
	kat = copy.deepcopy(tmpkat)

	def compute_gravity_tilt(tmpkat):
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
		#print " Wavefront tilt : %g nrad" % tilt
		return tilt

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
	kat.parseKatCode(code_WFS1)
	kat.parseKatCode(code_xaxis)
	kat.spo1.L=1000.0
	kat.ITM.ybeta=1.0e-10
	kat.ETM.ybeta=0.0
	t1=compute_gravity_tilt(kat)
	kat.ITM.ybeta=0.0
	kat.ETM.ybeta=-1.0e-10
	t2=compute_gravity_tilt(kat)

	kat = copy.deepcopy(tmpkat)
	kat.parseKatCode(code_WFS2)
	kat.parseKatCode(code_xaxis)
	kat.spo1.L=1.0e-9
	kat.ITM.ybeta=1.0e-10
	kat.ETM.ybeta=0.0
	t3=compute_gravity_tilt(kat)
	kat.ITM.ybeta=0.0
	kat.ETM.ybeta=-1.0e-10
	t4=compute_gravity_tilt(kat)
	print("	 WFS1 ITM {0:.4} nrad".format(t1*1e9))	  
	print("	 WFS2 ITM {0:.4} nrad".format(t3*1e9))	  
	print("	 WFS1 ETM {0:.4} nrad".format(t2*1e9))	  
	print("	 WFS2 ETM {0:.4} nrad".format(t4*1e9))	  
	return [t1,t2,t3,t4]
	
if __name__ == '__main__':
	main()

