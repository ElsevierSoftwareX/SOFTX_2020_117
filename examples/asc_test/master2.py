
from pykat import finesse
from pykat.commands import *
import pylab as pl
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
    
    
    out = tilt(kat)
    #(tilt_l, tilt_u) = asc_tilt.run(kat)
    
    kat.ETM.phi=result['phi_tuned']
    



def tilt(tmpkat):

    kat = copy.deepcopy(tmpkat)

    
    code_det = """
    beam PDrefl_car 0 nWFS2
    beam PDrefl_up 9M nWFS2
    beam PDrefl_low -9M nWFS2
    bp w0y y w0 nWFS2
    yaxis abs:deg
    """
    kat.parseKatCode(code_det)
    kat.noxaxis = True
    
    out = kat.run(printout=0,printerr=0)
    tilt_l = out.y[0]
    tilt_u = out.y[0]
    return (out)
    #return (tilt_l, tilt_u)

    
if __name__ == '__main__':
    main()

