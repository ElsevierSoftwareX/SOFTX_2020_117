from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import copy
from collections import namedtuple
from collections import OrderedDict
import pylab as pl
import shelve

import pykat
from pykat.components import *
from pykat.utilities.plotting.tools import printPDF
from pykat.external.progressbar import ProgressBar, ETA, Percentage, Bar
from pykat.utilities.plotting.tools import plot_setup
from pykat.optics.maps import *
from pykat.optics.gaussian_beams import HG_beam, beam_param
from pykat.optics.fft import *
from aligo import *

def main():
	print("""
	----------------------------------------
	""")


	# loading kat file to get parameters (if needed)
	global kat, out
	kat = pykat.finesse.kat()
	kat.loadKatFile('aligo_Xarm.kat')
	Lambda=kat.lambda0
	k = 2.0*np.pi/Lambda

	filename='fround-2014:12:22-14:17:37.npy'
	print(" --- loading data from file {0} ---".format(filename))
	global f_round
	f_round=np.load(filename)

	tmpresultfile = 'myshelf1.dat'
	# loading data saved by master.py
	try:
		tmpfile = shelve.open(tmpresultfile)
		result=tmpfile['result']
		tmpfile.close()
	except: raise Exception("Could not open temprary results file {0}".format(tmpresultfile))
        

	scan_start = 0.0
	scan_stop  = Lambda
	scan_points = 200
	global scan
	scan = np.linspace(scan_start, scan_stop, scan_points)

	# number of roundtrips
	global power
	N  = np.shape(f_round)[2]
	f_temp=np.zeros(np.shape(f_round[:,:,0]))
	power=np.zeros(scan_points,dtype=np.double)

	print(" --- performing cavity scan --- ")
	# This will take some time, let's show a progress bar
	p = ProgressBar(maxval=scan_points, widgets=["computing power:", Percentage(),"|", ETA(), Bar()])

	global phases, f_x, f_round
	ns=np.linspace(0.0, N-1, N)
	for i in range(scan_points):
		#f_temp[:,:]=0.0
		phases=np.exp(1j*2.0*k*scan[i]*ns)
		f_temp=np.sum(f_round*phases,axis=-1)
		#for n in range(N):
		#	f_temp = f_temp + np.multiply(f_round[:,:,n],np.exp(1j*k* 2.0*scan[i]*n));
		power[i] = field_power(f_temp,result['shape'])
		p.update(i)

	ax,fig=plot_setup()
	ax.plot(power)
	ax.set_yscale('log')
	pl.draw()
	pl.show(block=0)
	
def field_power(field, shape):
	return np.sum(np.abs(field)**2)*shape.xstep*shape.ystep;
	
if __name__ == '__main__':
    main()
