import copy
import sys
import scipy.optimize

def run(tmpkat):

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
        out = kat.run(printout=0,printerr=0)
        phase = (out.y[0]-out.y[1]-90)%360-180
        print '\r root finding: function value %g                    ' % phase ,
        sys.stdout.flush()
        return phase
    
    # do root finding
    xtol=1e-8
    (result, info)=scipy.optimize.bisect(carrier_resonance,0.0,40.0, xtol=xtol, maxiter=500, full_output=True)
    
    print ""
    if info.converged:
        print " Root has been found:"
        print " ETM phi %8f" % (result)
        print " (%d iterations, %g tolerance)" % (info.iterations, xtol)
        return result
    else:
        raise Exception(" Root has not been found")
        
        
if __name__=="__main__":
    run()
