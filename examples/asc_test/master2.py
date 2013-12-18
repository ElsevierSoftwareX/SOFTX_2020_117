
from pykat import finesse
from pykat.commands import *
import pylab as pl
import numpy as np
import shelve
import copy

def main():

    print """
    --------------------------------------------------------------
    Example file for using PyKat to automate Finesse simulations
    Finesse: http://www.gwoptics.org/finesse
    PyKat:   https://pypi.python.org/pypi/PyKat/
    
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
    global kat
    global out
    
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
    
    
    print "--------------------------------------------------------"
    print " 5. checking wavefront tilt for ITM/ETM tilt of 0.1nrad"

    #kat.PDrefl_p.enabled = False
    #kat.PDrefl_q.enabled = False

    kat.ETM.phi=result['phi_tuned']
    
    out = tilt(kat)
    #(tilt_l, tilt_u) = asc_tilt.run(kat)
    
    



def tilt(tmpkat):

    def compute_tilt(tmpkat):
        kat = copy.deepcopy(tmpkat)
        out = kat.run(printout=0,printerr=0)

        # compute x range in meters
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


    kat = copy.deepcopy(tmpkat)

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

    print " ETM ybeta 0.1nm"
    kat.ITM.ybeta=0.0
    kat.ETM.ybeta=1e-10
    (a3, a4) = compute_tilt(kat)


    
    print " WFS2:"
    print " ITM ybeta 0.1nm"
    kat = copy.deepcopy(tmpkat)
    kat.parseKatCode(code_comm)
    kat.parseKatCode(code_WFS2)
    kat.ITM.ybeta=1e-10
    kat.ETM.ybeta=0.0
    
    (a5, a6) = compute_tilt(kat)

    print " ETM ybeta 0.1nm"
    kat.ITM.ybeta=0.0
    kat.ETM.ybeta=1e-10
    (a6, a7) = compute_tilt(kat)

    out={}
    return (out)
    #return (tilt_l, tilt_u)

    
if __name__ == '__main__':
    main()

