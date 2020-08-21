# -*- coding: utf-8 -*-
# Author: Aaron Jones
# To test KatRun2D save and load
# Simple coupled cavity as test case

import numpy as np
import pykat
from pykat import finesse

filename = 'saveKatRun_output'
verbose = False
plot = False

def simulation_and_save(filename):
    kat = finesse.kat() # create a fresh cat object
    kat.verbose = verbose
    kat.parse("""
    l laser 1 0 n0
    s s_in 0.1 n0 n1
    m m_in 0.99 0.01 0 n1 n2
    s s_1 1 n2 n3
    m m_c 0.5 0.5 0 n3 n4
    s s_2 1 n4 n5
    m m_out 0.99 0.01 0 n5 n6
    pd trans n6
    xaxis m_in phi lin 0 180 360
    x2axis m_out phi lin 0 180 360
    """)
    out = kat.run()
    out.saveKatRun(filename+'.katrun')


def load(filename):
    out = finesse.KatRun2D.loadKatRun(filename+'.katrun')
    
    if verbose:
        print('x data:')
        print(out.x)
        print('x data type: \n\t',end='')
        print(type(out.x))
    
    if plot:
        from matplotlib import pyplot as plt 
        from matplotlib import colors
        
        fig, ax = plt.subplots()
        norm = colors.LogNorm(vmin=out['trans'].min(), 
                              vmax=out['trans'].max())
        pcm = ax.pcolor(out.x, out.y, out['trans'],norm=norm)
        cbar = fig.colorbar(pcm, ax=ax, extend='max')
        cbar.ax.set_ylabel('Transmitted Power [W]')
        ax.set_xlabel('ETM Tuning [deg]')
        ax.set_ylabel('PRM Tuning [deg]')
        fig.savefig(filename+'.png')
        plt.close(fig)
        print('Figure saved to: '+filename+'.png')

    assert isinstance(out.x,np.ndarray)
    assert isinstance(out.y,np.ndarray)
    assert isinstance(out['trans'],np.ndarray)

simulation_and_save(filename)
load(filename)