from itertools import combinations_with_replacement as combinations
from pykat.utilities.optics.gaussian_beams import beam_param, HG_beam
from pykat.exceptions import BasePyKatException

import maps
import os.path
import numpy as np
import pykat
import collections
import math

def makeCouplingMatrix(max_order):
    pass


def riemann_HG_knm(x, y, mode_in, mode_out, q1, q2, q1y=None, q2y=None, Axy=None):

    if Axy == None:
        Axy == np.ones((len(x), len(y)))
    
    if q1y == None:
        q1y = q1
        
    if q2y == None:
        q2y = q2
        
    if len(mode_in) != 2 or len(mode_out) != 2:
        raise BasePyKatException("Both mode in and out should be a container with modes [n m]")        
    
    Hg_in  = HG_beam(qx=q1, qy=q1y, n=mode_in[0], m=mode_in[1])
    Hg_out = HG_beam(qx=q2, qy=q2y, n=mode_out[0], m=mode_out[1])
    
    dx = abs(x[1] - x[0])
    dy = abs(y[1] - y[0])
    
    U1 = Hg_in.Unm(x,y)
    U2 = Hg_out.Unm(x,y).conjugate()
    
    return dx * dy * np.einsum('ij,ij', Axy, U1*U2)


def ROM_HG_knm(m, mode_in, mode_out, q1, q2, q1y=None, q2y=None):
    pass
    

def knmHG(couplings, surface_map, q1, q2, q1y=None, q2y=None, method="riemann"):
    
    couplings = np.array(couplings)
    
    a = couplings.size / 4.0
    
    if int(a) - a != 0:
        raise BasePyKatException("Iterator should be product of 4, each element of coupling array should be [n,m,n',m']")
        
    if "phase" in surface_map.type:
        Axy = np.exp(2j * 2 * math.pi / q1.wavelength * surface_map.data * surface_map.scaling)
    
    x = surface_map.x
    y = surface_map.y
    
    K = np.zeros((couplings.size/4,), dtype=np.complex128)
    
    it = np.nditer(couplings, flags=['refs_ok','f_index'])
    
    i = 0
    
    while not it.finished:
        try:
            
            mode_in = [int(it.next()), int(it.next())]
            mode_out = [int(it.next()), int(it.next())]
        
            if method == "riemann":
                K[i] = riemann_HG_knm(x, y, mode_in, mode_out, q1=q1, q2=q2, q1y=q1y, q2y=q2y, Axy=Axy)
            else:
                raise BasePyKatException("method value not accepted")
        
            i +=1
                 
        except StopIteration:
            break

    return K.reshape(couplings.shape[:-1])