from pykat import finesse
from pykat.commands import *
import pylab as pl
import scipy
#from scipy.optimize import minimize_scalar
import numpy as np
import shelve
import copy
import sys


def main():
    print """
    --------------------------------------------------------------
    Example file for using PyKat to automate Finesse simulations
    Finesse: http://www.gwoptics.org/finesse
    PyKat:   http://www.gwoptics.org/pykat
    
    The file runs through the various pykat files which are used
    to generate the Finesse results reported in the document:
    `Comparing Finesse simulations, analytical solutions and OSCAR 
    simulations of Fabry-Perot alignment signals', LIGO-T1300345
    
    This file is part of a collection.
    
    Andreas Freise 06.12.2013
    --------------------------------------------------------------
    """
    
    # shall we clear the workspace?
    # %reset -f

    # making these global during testing and debugging
    #global kat
    #global out
    
    kat = finesse.kat(tempdir=".",tempname="test")
    kat.verbose = False
    
    tmpresultfile = 'myshelf1.dat'
    
    # loading data saved by master.py
    kat.loadKatFile('asc_base2.kat')
    try:
        tmpfile = shelve.open(tmpresultfile)
        result=tmpfile['result']
        tmpfile.close()
    except: raise Exception("Could not open temprary results file {0}".format(tmpresultfile))
        
    # overwriting some variables
    kat.maxtem=3
    Lambda=1064.0e-9

    # disable PDH photo diode as we won't need it for most of this
    kat.PDrefl_p.enabled = False
    kat.PDrefl_q.enabled = False

    # simulating a tuned cavity
    kat.ETM.phi=result['phi_tuned']
    
    print "--------------------------------------------------------"
    print " 5. checking wavefronts for ITM/ETM tilt of 0.1nrad"
    tilt(kat)

    print "--------------------------------------------------------"
    print " 6. compute beam tilt from center of gravity calculation"
    gravity_tilt(kat)

    print "--------------------------------------------------------"
    print " 7. compute optimal demodulation phase of WFS1 and WFS2"

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
    kat.WFS1_I.phi[0]=WFS1_phase
    kat.WFS1_Q.phi[0]=WFS1_phase+90.0
    kat.WFS2_I.phi[0]=WFS2_phase
    kat.WFS2_Q.phi[0]=WFS2_phase+90.0
    result['WFS1_phase']=WFS1_phase
    result['WFS2_phase']=WFS2_phase

    print "--------------------------------------------------------"
    print " 8. compute ASC signal matrix at WFS1 and WFS2"
    signal = asc_signal(kat)
    
    print "--------------------------------------------------------"
    print " Saving results in temp. files to be read by master3.py"
    # re-enable PDH photo diode for savinf
    kat.PDrefl_p.enabled = True
    kat.PDrefl_q.enabled = True

    tmpkatfile = "asc_base3.kat"
    tmpresultfile = "myshelf2.dat"
    print " kat object saved in: {0}".format(tmpkatfile)
    print " current results saved in: {0}".format(tmpresultfile)
    # first the current kat file
    kat.saveScript(tmpkatfile)
    # now the result variables:
    tmpfile = shelve.open(tmpresultfile)
    tmpfile['result']=result
    tmpfile.close()


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
    out = kat.run(printout=0,printerr=0)
    WFS1_idx=out.ylabels.index("WFS1_I")
    WFS2_idx=out.ylabels.index("WFS2_I")
    signal[0,0] = out.y[WFS1_idx]
    signal[1,0] = out.y[WFS2_idx]

    kat.ITM.ybeta=0.0
    kat.ETM.ybeta=-1e-10
    out = kat.run(printout=0,printerr=0)
    signal[0,1] = out.y[WFS1_idx]
    signal[1,1] = out.y[WFS2_idx]
    signal = signal *1e10
    sensors=('WFS1', 'WFS2')
    mirrors=('ITM', 'ETM')
    print "  ASC Matrix:"
    for i in range(2):
        print "  ", sensors[i], " ",
        for j in range(2):
            print "%12.10g" % signal[i,j],
        print mirrors[i]
    return signal
    
    
def asc_phases(tmpkat):
    kat = copy.deepcopy(tmpkat)
    
    kat.parseKatCode('yaxis abs')
    kat.noxaxis = True
    kat.maxtem=1

    def demod_phase1(x):
        kat.WFS1_I.phi[0]=x[0]
        out = kat.run(printout=0,printerr=0)
        WFS1_idx=out.ylabels.index("WFS1_I")
        signal = out.y[WFS1_idx]
        print '\r minimising: function value %g                    ' % signal ,
        sys.stdout.flush()
        return -1*abs(signal)

    def demod_phase2(x):
        kat.WFS2_I.phi[0]=x[0]
        out = kat.run(printout=0,printerr=0)
        WFS2_idx=out.ylabels.index("WFS2_I")
        signal = out.y[WFS2_idx]
        print '\r minimising: function value %g                    ' % signal ,
        sys.stdout.flush()
        return -1*abs(signal)

    kat.ITM.ybeta=1e-10
    kat.ETM.ybeta=0.0
    # minimize_scaler is only available in newer scipy versions
    #res = minimize_scalar(demod_phase1, method='brent')
    if int(scipy.version.version.split('.')[1])<11:
        from scipy.optimize import fmin
        fmin(demod_phase1, x0, xtol=1e-8)
    else:
        from scipy.optimize import minimize
        res = minimize(demod_phase1, 0, method='nelder-mead', options={'xtol':1e-8,'disp': False})
    WFS1_phase = res.x[0]
    print ""
    print " WFS1 demod phase : %.10g deg" % WFS1_phase
     
    kat.ITM.ybeta=0.0
    kat.ETM.ybeta=-1e-10
    # minimize_scaler is only available in newer scipy versions
    #res = minimize_scalar(demod_phase2, method='brent')
    if int(scipy.version.version.split('.')[1])<11:
        from scipy.optimize import fmin
        fmin(demod_phase2, x0, xtol=1e-8)
    else:
        from scipy.optimize import minimize
        res = minimize(demod_phase2, 0, method='nelder-mead', options={'xtol':1e-8,'disp': False})
    WFS2_phase = res.x[0]
    print ""
    print " WFS2 demod phase : %.10g deg" % WFS2_phase
    return(WFS1_phase, WFS2_phase)    
    
def gravity_tilt(tmpkat):
    kat = copy.deepcopy(tmpkat)

    def compute_gravity_tilt(tmpkat):
        kat = copy.deepcopy(tmpkat)
        out = kat.run(printout=0,printerr=0)

        y1 = out["b1"]
        y2 = out["b1_1k"]
        # shift of beam center  on detector 1 (as m/w0y)
        x1 = np.sum(out.x*y1)/np.sum(y1) 
        # shift of beam center  on detector 2 (as m/w0y)
        x2 = np.sum(out.x*y2)/np.sum(y2)
        # calibrate this in meter by mutliplying with w0y
        # and compute the angle geometrically        
        w0=out["w0y"][0]
        detector_distance = 1000.0
        tilt=w0*(x2-x1)/detector_distance
        print " Wavefront tilt : %g nrad" % tilt

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
    print " WFS1:"
    print " ITM ybeta 0.1nm"
    kat.parseKatCode(code_WFS1)
    kat.parseKatCode(code_xaxis)
    kat.spo1.L=1000.0
    kat.ITM.ybeta=1e-10
    kat.ETM.ybeta=0.0
    compute_gravity_tilt(kat)
    print " ETM ybeta -0.1nm"
    kat.ITM.ybeta=0.0
    kat.ETM.ybeta=-1e-10
    compute_gravity_tilt(kat)

    print " WFS1:"
    print " ITM ybeta 0.1nm"
    kat = copy.deepcopy(tmpkat)
    kat.parseKatCode(code_WFS2)
    kat.parseKatCode(code_xaxis)
    kat.spo1.L=1.0e-9
    kat.ITM.ybeta=1e-10
    kat.ETM.ybeta=0.0
    compute_gravity_tilt(kat)
    print " ETM ybeta -0.1nm"
    kat.ITM.ybeta=0.0
    kat.ETM.ybeta=-1e-10
    compute_gravity_tilt(kat)

    
def tilt(tmpkat):
    kat = copy.deepcopy(tmpkat)
    
    def compute_tilt(tmpkat):
        kat = copy.deepcopy(tmpkat)
        out = kat.run(printout=0,printerr=0)

        # compute data x range in meters
        beamsize = out["w0y"][0,0] 
        xrange = beamsize*(out.x.max()-out.x.min())
        stepsize=xrange/(len(out.x)-1)
        print " Beamsize %e m" % beamsize
        print " Measurement range: %e m, stepszie: %e m" % (xrange, stepsize)
        # compute difference in angle between wavefront of carrier and sidebands
        diff_l = (out["PDrefl_low"][:,1]-out["PDrefl_car"][:,1])/stepsize
        diff_u = (out["PDrefl_up"][:,1]-out["PDrefl_car"][:,1])/stepsize
        tilt_l = diff_l[1:-1]-diff_l[0:-2]
        tilt_u = diff_u[1:-1]-diff_u[0:-2]
        print " Tilt (upper  - car), mean: %e m/deg, stddev %e m/deg" % (np.mean(tilt_u), np.std(tilt_u))
        print " Tilt (lower  - car), mean: %e m/deg, stddev %e m/deg" % (np.mean(tilt_l), np.std(tilt_l))
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

    print " WFS1:"
    print " ITM ybeta 0.1nm"
    kat.parseKatCode(code_comm)
    kat.parseKatCode(code_WFS1)
    kat.ITM.ybeta=1e-10
    kat.ETM.ybeta=0.0
    (a1, a2) = compute_tilt(kat)

    print " ETM ybeta -0.1nm"
    kat.ITM.ybeta=0.0
    kat.ETM.ybeta=-1e-10
    (a3, a4) = compute_tilt(kat)
    
    print " WFS2:"
    print " ITM ybeta 0.1nm"
    kat = copy.deepcopy(tmpkat)
    kat.parseKatCode(code_comm)
    kat.parseKatCode(code_WFS2)
    kat.ITM.ybeta=1e-10
    kat.ETM.ybeta=0.0
    (a5, a6) = compute_tilt(kat)

    print " ETM ybeta -0.1nm"
    kat.ITM.ybeta=0.0
    kat.ETM.ybeta=-1e-10
    (a6, a7) = compute_tilt(kat)

    return 
    
if __name__ == '__main__':
    main()

