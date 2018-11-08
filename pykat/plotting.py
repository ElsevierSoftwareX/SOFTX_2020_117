# -*- coding: utf-8 -*-
"""
Created on Sat Feb 02 10:35:04 2013

@author: Daniel
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

__DPI__ = None


def in_ipython():
    try:
        __IPYTHON__
    except NameError:
        return False
    else:
        return True


# `mode` should be either display (showing in windows)
# or paper (saving for paper/report/etc.)
def init_pykat_plotting(mode="display", dpi=None, fmts=None):
    import matplotlib as mpl
    import pykat.style

    if fmts is None:
        fmts = ['svg']

    if in_ipython():
        try:
            from IPython.display import set_matplotlib_formats
            ipy = get_ipython()
            ipy.magic("matplotlib inline")
            set_matplotlib_formats(*fmts)
        except NameError:
            pass

    if mode == "display":
        pykat.style.use(["default"])
    elif mode == "paper":
        mpl.use("pgf")
        pykat.style.use(["default", "paper"])
    else:
        raise (BaseException(
            "Plotting mode must be either 'display' or 'paper'."))

    if dpi is None:
        __DPI__ = mpl.rcParams['figure.dpi']
    else:
        __DPI__ = int(dpi)
        mpl.rcParams.update({'figure.dpi': __DPI__})
        mpl.rcParams.update({'savefig.dpi': __DPI__})

    if (mpl.__version__ < '1.5'):
        # prop_cycle doesn't exist pre-1.5, and unknown rcParams are ignored,
        # so we have to hardcode the color cycle
        mpl.rcParams.update({'axes.color_cycle': ['0000FF', 'FF0000', '000000', '00FF00', 'FF00FF']})


def figure(width=None, height=None, textwidth=None, **kwargs):
    """
    Options:
        width: 'full', 'half' (0.49*textwidth) (default: None (use current
        value))
        height: relative height to width (default: None (use current value))
        textwidth: Width of text in inches (default: None (use current value))
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    fig_size = mpl.rcParams["figure.figsize"]

    if width is None:
        width = "full"
    if height is None:
        height = fig_size[1] / fig_size[0]
    if textwidth is None:
        textwidth = fig_size[0]

    if width == "full":
        fig_size[0] = textwidth
    elif width == "half":
        fig_size[0] = 0.49 * textwidth
        fig_size[1] = 0.49 * textwidth * height
    else:
        raise (BaseException("width must be either 'full', 'half' or None."))

    fig = plt.figure(figsize=fig_size, dpi=__DPI__)

    return fig


def subplot(*args, **kwargs):
    import matplotlib.pyplot as plt

    ax = plt.subplot(*args, **kwargs)

    return ax
