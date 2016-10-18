from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import matplotlib.pyplot as plt

from pykat.optics.gaussian_beams import HG_mode, beam_param
from pykat.optics.fft import *
import numpy as np
import scipy

def main():
    print("""
    --------------------------------------------------------------
    Example file for using PyKat http://www.gwoptics.org/pykat

    Generate a HG11 mode from a HG00 mode with a simple
    phase plate.
    
    Andreas Freise, 18.10.2016    
    --------------------------------------------------------------
    """)
    plt.close('all')
    # wavelength
    Lambda = 1064.0E-9 
    # distance to propagate/focal length of lens
    D = 4

    ######## Generate Grid stucture required for FFT propagation ####
    xpoints = 512
    ypoints = 512
    xsize = 0.05
    ysize = 0.05
    # Apply offset such that the center of the beam lies in the
    # center of a grid tile
    xoffset = -0.5*xsize/xpoints
    yoffset = -0.5*ysize/ypoints

    shape = grid(xpoints, ypoints, xsize, ysize, xoffset, yoffset)
    x = shape.xaxis
    y = shape.yaxis

    ######## Generates input beam ################
    gx = beam_param(w0=2e-3, z=0)
    gy = beam_param(w0=2e-3, z=0)
    beam = HG_mode(gx,gy,0,0)
    global field, laser
    field = beam.Unm(x,y) 
 
    ####### Apply phase plate #######################################

    plate = np.ones([xpoints, ypoints])
    plate[0:xpoints//2, 0:ypoints//2]=-1.0
    plate[xpoints//2:xpoints, ypoints//2:ypoints]=-1.0
    
    field2 = field*plate;

    ####### Propagates the field by FFT ##############################
    field2 = FFT_propagate(field2,shape,Lambda,D,1) 


    # maybe apply a thin lens
    f=D
    field2 = apply_thin_lens(field2, shape, Lambda, f) 
    field2 = FFT_propagate(field2,shape,Lambda,D,1)
    
    midx=(xpoints)//2
    midy=(ypoints)//2
    off1=50
    off2=50

    # plot hand tuned for certain ranges and sizes, not automtically scaled
    fig=plt.figure(110)
    fig.clear()
    plt.subplot(1, 3, 1)
    plt.imshow(abs(field))
    plt.xlim(midx-off1,midx+off1)
    plt.ylim(midy-off1,midy+off1)
    plt.draw()
    plt.subplot(1, 3, 2)
    plt.imshow(plate)
    #pl.xlim(midx-off2,midx+off2)
    #pl.ylim(midy-off2,midy+off2)
    plt.draw()
    plt.subplot(1, 3, 3)
    plt.imshow(abs(field2))
    plt.xlim(midx-off2,midx+off2)
    plt.ylim(midy-off2,midy+off2)
    plt.draw()
    if in_ipython():
        plt.show(block=0)
    else:
        plt.show(block=1)


# testing if the script is run from within ipython
def in_ipython():
    try:
        __IPYTHON__
    except NameError:
        return False
    else:
        return True

if __name__ == '__main__':
    main()

