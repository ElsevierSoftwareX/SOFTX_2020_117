# -*- coding: utf-8 -*-
"""
Created on Sat Feb 02 10:35:04 2013

@author: Daniel
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import matplotlib

BACKEND = 'Qt4Agg'
matplotlib.use(BACKEND)

from matplotlib import rc
import matplotlib.pyplot as plt

mainpid = -1

def plot1D(run, title=""):
    
    rc('font', **pp.font)
    rc('xtick',labelsize=pp.TICK_SIZE)
    rc('ytick',labelsize=pp.TICK_SIZE)
    rc('text', usetex=pp.USETEX)
    rc('axes', labelsize = pp.LABEL_SIZE)
    
    fig=plt.figure()
    fig.set_size_inches(pp.fig_size)
    fig.set_dpi(pp.FIG_DPI)
        
    ax1 = fig.add_subplot(111)
    ax1.set_xlim(np.min(run.x),np.max(run.x))
    traces = ax1.plot(run.x,run.y)
    ax1.grid(pp.GRID)
    
    ax1.set_xlabel(run.xlabel)
    legends = run.ylabels
    ax1.legend(traces, legends, loc=0, shadow=pp.SHADOW,prop={'size':pp.LEGEND_SIZE})

    if pp.PRINT_TITLE:
        plt.title(title)
        
    if pp.SCREEN_TITLE:
        fig.canvas.manager.set_window_title(title)
    else:
        fig.canvas.manager.set_window_title('')
        
    #plt.ion()
    plt.show()
    
class pp():
    # set some gobal settings first
    BACKEND = 'Qt4Agg' # matplotlib backend
    FIG_DPI=90 # DPI of on sceen plot
    # Some help in calculating good figure size for Latex
    # documents. Starting with plot size in pt,
    # get this from LaTeX using \showthe\columnwidth
    fig_width_pt = 484.0
    inches_per_pt = 1.0/72.27  # Convert TeX pt to inches
    golden_mean = (np.sqrt(5)-1.0)/2.0   # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size = [fig_width,fig_height]
    # some plot options:
    LINEWIDTH = 1 # linewidths of traces in plot
    AA = True # antialiasing of traces
    USETEX = False # use Latex encoding in text
    SHADOW = False # shadow of legend box
    GRID = True # grid on or off
    # font sizes for normal text, tick labels and legend
    FONT_SIZE = 10 # size of normal text
    TICK_SIZE = 10 # size of tick labels
    LABEL_SIZE = 10 # size of axes labels
    LEGEND_SIZE = 10 # size of legend
    # font family and type
    font = {'family':'sans-serif','sans-serif':['Helvetica'],'size':FONT_SIZE}
    DPI=300 # DPI for saving via savefig
    # print options given to savefig command:
    print_options = {'dpi':DPI, 'transparent':True, 'bbox_inches':'tight', 'pad_inches':0.1}
    # for Palatino and other serif fonts use:
    #font = {'family':'serif','serif':['Palatino']}
    SCREEN_TITLE = True # show title on screen?
    PRINT_TITLE = False # show title in saved file?
