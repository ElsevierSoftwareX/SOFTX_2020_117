import copy
import sys
import scipy.optimize

from pykat import finesse

def run(tmpkat):

    kat = copy.deepcopy(tmpkat)
    
    code_det = """
    pd1 PDrefl_q 9M 90 nWFS1
    %scale 2 PDrefl_q
    """
    
    kat.parseKatCode(code_det)
    kat.noxaxis= True

    # function for root finding
    def PD_q_test(x):
        kat.PDrefl_q.phase[0]=x
        out = kat.run(printout=0,printerr=0)
        print '\r root finding: function value %g                    ' % out.y,
        sys.stdout.flush()
        return out.y

    # do root finding
    xtol=1e-8
    (result, info)=scipy.optimize.bisect(PD_q_test,80.0,100.0, xtol=xtol, maxiter=500, full_output=True)

    print ""
    if info.converged:
        p_phase=result-90.0
        q_phase=result
        print " Root has been found:"
        print " p_phase %8f" % (p_phase)
        print " q_phase %8f" % (q_phase)
        print " (%d iterations, %g tolerance)" % (info.iterations, xtol)
        return (p_phase, q_phase)
    else:
        raise Exception("Root has not been found")
        

    
