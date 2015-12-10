# -*- coding: utf-8 -*-
"""
Created on Sat Feb 02 10:35:04 2013

@author: Daniel
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# Should be either display (showing in windows) or paper (saving for paper/report/etc.)
__mode__ = None
__DPI__ = None

def in_ipython():
    try:
        cfg = get_ipython() 
        return True
    except NameError:
        return False
        
        
def init_pykat_plotting(mode="display", dpi=100):
    import matplotlib as mpl
    
    __DPI__ = int(dpi)
    
    if in_ipython():
        try:
            from IPython.display import set_matplotlib_formats
            set_matplotlib_formats('pdf', 'svg')
            ipy = get_ipython()
            ipy.magic("matplotlib inline")
        except:
            pass
    
    if mode == "display":
        __mode__ = mode
        
    elif mode == "paper":
        __mode__ = mode
        
        mpl.use("pgf")

        pgf_with_pdflatex = {
            "pgf.texsystem": "pdflatex",
            "pgf.preamble": [
                 r"\\usepackage{amsmath, amssymb}",
                 r"\\usepackage{mathtools, siunitx}" ,
                 r"\\usepackage{amsmath}",
                 r"\\usepackage[utf8x]{inputenc}",
                 r"\\usepackage[T1]{fontenc}"
                 ]
        }

        mpl.rcParams.update(pgf_with_pdflatex)
    else:
        raise(BaseException("Plotting mode must be either 'display' or 'paper'."))

    if (mpl.__version__ < '1.5'):
        mpl.rcParams['axes.color_cycle'] = ['b', 'r', 'k', 'g', 'c', 'm', 'y']
    else:
        mpl.rcParams['axes.prop_cycle']=mpl.cycler('color', ['b', 'r', 'k', 'g', 'c', 'm', 'y'])
    mpl.rcParams['lines.linewidth'] = 1.2
    mpl.rcParams.update({"figure.figsize": (6, 3.708)})
    mpl.rcParams.update({'font.size': 11})
    mpl.rcParams.update({'figure.dpi': __DPI__})
    mpl.rcParams.update({'savefig.dpi': __DPI__})
    mpl.rcParams.update({'font.family': "serif"})
    mpl.rcParams.update({'axes.grid': True})
    mpl.rcParams.update({'axes.axisbelow': True})
    mpl.rcParams.update({'grid.linewidth': 0.25})
    mpl.rcParams.update({'grid.linestyle': ":"})
    mpl.rcParams.update({'grid.color': (0.6,0.6,0.6,1)})
    mpl.rcParams.update({'savefig.bbox': "tight"})
    mpl.rcParams.update({'savefig.pad_inches': 0.05})
    mpl.rcParams.update({'xtick.labelsize': "small"})
    mpl.rcParams.update({'ytick.labelsize': "small"})
    mpl.rcParams.update({'axes.formatter.useoffset': False})

def figure(width="full", height=0.618, textwidth=6, **kwargs):
    """
    Options:
        width: 'full', 'half' (0.49*textwidth) (default: full)
        height: relative height to width (default: 1/golden ratio = 0.618)
        textwidth: Width of text in inches (default: 9.84252 in = 25 cm )
    """
    import matplotlib.pyplot as plt
    
    if width == "full":
        fig_size = [textwidth, textwidth*height]
    elif width == "half":
        fig_size = [textwidth*0.49, textwidth*height*0.49]
    else:
        raise(BaseException("width must be either 'full' or 'half'."))
        
    fig = plt.figure(figsize=fig_size, dpi=__DPI__)
    
    return fig
    
    
    
    
