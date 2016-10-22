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
#from pykat.tools.plotting import printPDF
from pykat.external.progressbar import ProgressBar, ETA, Percentage, Bar, Timer
import pykat.plotting 
from pykat.optics.maps import *
from pykat.optics.gaussian_beams import HG_mode, beam_param
from pykat.optics.fft import *
from aligo import *

def main():
    print("""
    ----------------------------------------
    """)


    # loading kat file to get parameters (if needed)
    global kat, out
    kat = pykat.finesse.kat()
    kat.verbose = False
    kat.loadKatFile('aligo_Xarm.kat')
    Lambda=kat.lambda0
    k = 2.0*np.pi/Lambda

    #filename='fround_mode_matched_no_map.npy'
    #filename='fround_mode_matched_10map.npy'
    filename='fround-2014:12:29-19:55:28.npy'
    filename='fround-2016:10:21-16:10:14.npy'
    print(" --- loading data from file {0} ---".format(filename))
    global f_round
    f_round=np.load(filename)

    tmpresultfile = 'myshelf1.dat'
    # loading additional data saved by previous file
    try:
        tmpfile = shelve.open(tmpresultfile)
        result=tmpfile['result']
        tmpfile.close()
    except: raise Exception("Could not open temprary results file {0}".format(tmpresultfile))
        
    scan_start = 0.0
    scan_stop  = Lambda
    scan_points = 100
    global scan
    scan = np.linspace(scan_start, scan_stop, scan_points)

    # number of roundtrips
    global power
    N  = np.shape(f_round)[2]
    f_temp=np.zeros(np.shape(f_round[:,:,0]))
    power=np.zeros(scan_points,dtype=np.double)

    print(" --- performing cavity scan --- ")
    # This will take some time, let's show a progress bar
    p = ProgressBar(maxval=scan_points, widgets=["power:", Percentage(),"|", Timer(), "|", ETA(), Bar()])

    global phases
    ns=np.linspace(0.0, N-1, N)
    for i in range(scan_points):
        phases=np.exp(1j*2.0*k*scan[i]*ns)
        f_temp=np.sum(f_round*phases,axis=-1)
        power[i] = field_power(f_temp,result['shape'])
        p.update(i)
    p.finish()

    fileName, fileExtension = os.path.splitext(filename)
    txtfile='power_{0}.txt'.format(fileName)
    np.savetxt(txtfile, power, fmt='%.18e', delimiter=' ')

    # plot scan 
    #ax,fig=plot_setup()
    fig=pykat.plotting.figure()
    ax=fig.add_subplot(111)
    ax.plot(power)
    ax.set_yscale('log')
    pl.draw()
    pl.show(block=0)
    
def field_power(field, shape):
    return np.sum(np.abs(field)**2)*shape.xstep*shape.ystep;
    
if __name__ == '__main__':
    main()
