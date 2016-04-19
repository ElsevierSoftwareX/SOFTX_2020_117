from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from pykat import finesse
from pykat.commands import *
import numpy as np
import pickle
import copy
import sys
import shutil

def main():

    print("""
    --------------------------------------------------------------
    Example file for using PyKat to automate Finesse simulations
    Finesse: http://www.gwoptics.org/finesse
    PyKat:   http://www.gwoptics.org/pykat

    The file runs through the various Finesse simulations
    to generate the Finesse results reported in the document:
    `Comparing Finesse simulations, analytical solutions and OSCAR 
    simulations of Fabry-Perot alignment signals', LIGO-T1300345,
    freely available online: http://arxiv.org/abs/1401.5727

    Run this file after master2.py to create data which can be
    plotted using master3_plot.py. Results are saved after 
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

    kat.PDrefl_q.enabled = False
    kat.WFS1_Q.enabled = False
    kat.WFS2_Q.enabled = False

    print("--------------------------------------------------------")
    print(" 9. ASC signals for large misalignments (ITM)")
    asc_large(kat, 'ITM')

    print("--------------------------------------------------------")
    print(" 10. ASC signals for large misalignments (ETM)")
    asc_large(kat, 'ETM')


#-----------------------------------------------------------------------------------

def asc_large(tmpkat, mir_name):
    kat = copy.deepcopy(tmpkat)

    assert(mir_name == 'ITM' or mir_name == 'ETM')
    
    code_lock = """
    set err PDrefl_p re
    lock z $err 900 1p
    put* ETM phi $z
    noplot z
    """
        
    kat.parseKatCode(code_lock)
    kat.parseKatCode('yaxis abs')
    kat.parseKatCode('xaxis {0} ybeta lin 0 1u 100'.format(mir_name))
    maxtems = [1, 3, 5]
    #kat.verbose=1
    xscale = 1e6
    yscale = 1e6
    #global out
    tmpfilename = "datashelf_{0}.dat".format(mir_name)
    backupname = "datashelf_{0}.dat.bck".format(mir_name)
    out={}
    done_maxtems = []
    
    for tem in maxtems:
        done_maxtems.append(tem)
        print(" Calculating maxtem = %d " % tem)
        kat.maxtem = tem
        out[str(tem)] = kat.run()
        import os.path
        if os.path.isfile(tmpfilename):
            shutil.copyfile(tmpfilename, backupname)
        print(" current results saved in: {0}".format(tmpfilename))
        with open(tmpfilename, 'wb') as handle:
            pickle.dump({ "out": out, "maxtems": done_maxtems}, handle)    
    
if __name__ == '__main__':
    main()
