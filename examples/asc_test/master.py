from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from pykat import finesse
#from pykat.commands import *
import copy
import pickle
import sys
import scipy.optimize


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
    shown the document's sections 3 and 4 and saves temporary
    data and a new Finesse input file to be read by master2.py.
    
    Andreas Freise 16.01.2014
    --------------------------------------------------------------
    """)   
    
    # for debugging we might need to see the temporay file:
    global kat
    kat = finesse.kat(tempdir=".",tempname="test")
    kat.verbose = False
    kat.loadKatFile('asc_base.kat')
    kat.maxtem=3
    Lambda=1064.0e-9
    result = {}
    # defining variables as global for debugging
    #global kat, out, result
    
    print("--------------------------------------------------------")
    print(" 1. tunes ETM position to find resonance")
    kat.ETM.phi=resonance(kat)
    
    print("--------------------------------------------------------")
    print(" 2. print sideband and carrier powers/amplitudes")
    powers(kat)
    
    print("--------------------------------------------------------")
    print(" 3. determine the optimal phase for the PDH signal")
    (result['p_phase'], result['q_phase']) = pd_phase(kat)
    
    # setting demodulation phase
    code_det = """
    pd1 PDrefl_p 9M 0 nWFS1
    scale 2 PDrefl_p
    pd1 PDrefl_q 9M 90 nWFS1
    scale 2 PDrefl_q
    """
    kat.parseKatCode(code_det)
    kat.PDrefl_p.phase1=result['p_phase']
    kat.PDrefl_q.phase1=result['q_phase']
    
    print("--------------------------------------------------------")
    print(" 4. adding a 0.1nm offset to ETM and compute PDH signal")
    result['phi_tuned']=float(kat.ETM.phi)
    result['phi_detuned'] = result['phi_tuned'] + 0.1*360.0/1064.0
    
    kat.ETM.phi=result['phi_detuned']
    print(" new ETM phi tuning = %g " % kat.ETM.phi)

    (result['pd_p'], result['pd_q']) = pd_signal(kat)
    print(" PDH inphase     = %e " % result['pd_p'])
    print(" PDH quadrtature = %e " % result['pd_q'])
    
    print("--------------------------------------------------------")
    print(" Saving results in temp. files to be read by master2.py")
    tmpkatfile = "asc_base2.kat"
    tmpresultfile = "myshelf1.dat"
    print(" kat object saved in: {0}".format(tmpkatfile))
    print(" current results saved in: {0}".format(tmpresultfile))
    # first the current kat file
    kat.saveScript(tmpkatfile)
    with open(tmpresultfile, 'wb') as handle:
        pickle.dump(result, handle)
    
#---------------------------------------------------------------------------

def pd_signal(tmpkat):

    kat = copy.deepcopy(tmpkat)
    code1="""
        pd cav nITM2
        yaxis abs
        """
    kat.parseKatCode(code1)
    kat.noxaxis = True
    global out
    out = kat.run()
    print(" Cavity power: {0:.6f}W".format(float(out['cav'])))
    return (float(out['PDrefl_p']), float(out['PDrefl_q']))
    
def pd_phase(tmpkat):

    kat = copy.deepcopy(tmpkat)
    code_det = """
    pd1 PDrefl_q 9M 90 nWFS1
    """
    
    kat.parseKatCode(code_det)
    kat.noxaxis= True

    # function for root finding
    def PD_q_test(x):
        kat.PDrefl_q.phase1=x
        out = kat.run()
        print('\r root finding: function value {0:<16g}'.format(float(out.y)), end='')
        sys.stdout.flush()
        return float(out.y)

    # do root finding
    xtol=1e-8

    (result, info)=scipy.optimize.bisect(PD_q_test,80.0,100.0, xtol=xtol, maxiter=500, full_output=True)

    print("")
    if info.converged:
        p_phase=result-90.0
        q_phase=result
        print(" Root has been found:")
        print(" p_phase %8f" % (p_phase))
        print(" q_phase %8f" % (q_phase))
        print(" (%d iterations, %g tolerance)" % (info.iterations, xtol))
        return (p_phase, q_phase)
    else:
        raise Exception("Root has not been found")
        

def powers(tmpkat):

    kat = copy.deepcopy(tmpkat)
    
    code1 = """
    ad EOM_up 9M nEOM1
    ad EOM_low -9M nEOM1
    pd cav_pow nITM2
    ad cav_c 0 nITM2
    ad WFS1_u  9M nWFS1
    ad WFS1_l -9M nWFS1
    ad WFS1_c  0  nWFS1
    ad WFS2_u  9M nWFS2
    ad WFS2_l -9M nWFS2
    ad WFS2_c	0 nWFS2
    noxaxis
    """

    kat.parseKatCode(code1)

    global out
    out = kat.run()
    for i in range(len(out.y[0])):
        print(" %8s: %.4e" % (out.ylabels[i], out.y[0,i]))


def resonance(tmpkat):
    kat = copy.deepcopy(tmpkat)
    
    code1 = """
    ad carr2 0 nITM1*
    ad carr3 0 nITM2
    yaxis deg
    """
    kat.parseKatCode(code1)
    kat.noxaxis = True
    
    # function for root finding
    def carrier_resonance(x):
        kat.ETM.phi=x
        out = kat.run()
        phase = (out.y[0,0]-out.y[0,1]-90)%360-180.0
        print('\r root finding: function value {0:<16g}'.format(float(phase)), end='')
        sys.stdout.flush()
        return phase
    
    # do root finding
    xtol=1e-8
    (result, info)=scipy.optimize.bisect(carrier_resonance,0.0,40.0, xtol=xtol, maxiter=500, full_output=True)
    
    print("")
    if info.converged:
        print(" Root has been found:")
        print(" ETM phi %8f" % (result))
        print(" (%d iterations, %g tolerance)" % (info.iterations, xtol))
        return result
    else:
        raise Exception(" Root has not been found")
        

if __name__ == '__main__':
    main()

