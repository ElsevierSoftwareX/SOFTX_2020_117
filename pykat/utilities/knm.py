from itertools import combinations_with_replacement as combinations
from pykat.utilities.optics.gaussian_beams import beam_param, HG_beam
from pykat.exceptions import BasePyKatException

import maps
import os.path
import numpy as np
import pykat
import collections
import math
from romhom import u_star_u

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


def __gen_ROM_HG_knm_cache(weights, couplings, q1, q2, q1y=None, q2y=None):

    if q1y == None:
        q1y = q1
        
    if q2y == None:
        q2y = q2
        
    it = np.nditer(couplings, flags=['refs_ok','f_index'])
    
    cache = {}
    
    cache["w_ij_Q1Q3"] = weights.w_ij_Q1 + weights.w_ij_Q3
    cache["w_ij_Q2Q4"] = weights.w_ij_Q2 + weights.w_ij_Q4
    cache["w_ij_Q1Q2"] = weights.w_ij_Q1 + weights.w_ij_Q2
    cache["w_ij_Q1Q4"] = weights.w_ij_Q1 + weights.w_ij_Q4
    cache["w_ij_Q2Q3"] = weights.w_ij_Q2 + weights.w_ij_Q3
    cache["w_ij_Q3Q4"] = weights.w_ij_Q3 + weights.w_ij_Q4
    cache["w_ij_Q1Q2Q3Q4"] = weights.w_ij_Q1 + weights.w_ij_Q3 + weights.w_ij_Q2 + weights.w_ij_Q4
    
    while not it.finished:
        try:
            mode_in = [int(it.next()), int(it.next())]
            mode_out = [int(it.next()), int(it.next())]
            
            strx = "x[%i,%i]" % (mode_in[0], mode_out[0])
            stry = "y[%i,%i]" % (mode_in[1], mode_out[1])
            
            if strx not in cache:
                cache[strx] = u_star_u(q1.z,   q2.z,  q1.w0,  q2.w0, mode_in[0], mode_out[0], weights.nodes)   
            
            if q1 == q1y and q2 == q2y:
                cache[stry] = cache[strx]
            elif stry not in cache:
                cache[stry] = u_star_u(q1y.z, q2y.z, q1y.w0, q2y.w0, mode_in[1], mode_out[1], weights.nodes)
            
            
            
        except StopIteration:
            break
    
    return cache
    
def ROM_HG_knm(weights, mode_in, mode_out, q1, q2, q1y=None, q2y=None, cache=None):
    if q1y == None:
        q1y = q1

    if q2y == None:
        q2y = q2
    
    # x modes
    n = mode_in[0]
    m = mode_out[0]

    # y modes
    npr = mode_in[1]
    mpr = mode_out[1]

    if cache == None:
        u_x_nodes = u_star_u(q1.z,   q2.z,  q1.w0,  q2.w0, n,   m,   weights.nodes)   
        u_y_nodes = u_star_u(q1y.z, q2y.z, q1y.w0, q2y.w0, npr, mpr, weights.nodes)
        
        w_ij_Q1Q3 = weights.w_ij_Q1 + weights.w_ij_Q3
        w_ij_Q2Q4 = weights.w_ij_Q2 + weights.w_ij_Q4
        w_ij_Q1Q2 = weights.w_ij_Q1 + weights.w_ij_Q2
        w_ij_Q1Q4 = weights.w_ij_Q1 + weights.w_ij_Q4
        w_ij_Q2Q3 = weights.w_ij_Q2 + weights.w_ij_Q3
        w_ij_Q3Q4 = weights.w_ij_Q3 + weights.w_ij_Q4
        w_ij_Q1Q2Q3Q4 = weights.w_ij_Q1 + weights.w_ij_Q3 + weights.w_ij_Q2 + weights.w_ij_Q4
    else:
        strx = "x[%i,%i]" % (mode_in[0], mode_out[0])
        stry = "y[%i,%i]" % (mode_in[1], mode_out[1])
        
        u_x_nodes = cache[strx]
        u_y_nodes = cache[stry]
        
        w_ij_Q1Q3 = cache["w_ij_Q1Q3"]
        w_ij_Q2Q4 = cache["w_ij_Q2Q4"]
        w_ij_Q1Q2 = cache["w_ij_Q1Q2"]
        w_ij_Q1Q4 = cache["w_ij_Q1Q4"]
        w_ij_Q2Q3 = cache["w_ij_Q2Q3"]
        w_ij_Q3Q4 = cache["w_ij_Q3Q4"]
        w_ij_Q1Q2Q3Q4 = cache["w_ij_Q1Q2Q3Q4"]
    
    u_xy_nodes = np.outer(u_x_nodes, u_y_nodes)
    
    n_mod_2 = n % 2
    m_mod_2 = m % 2
    npr_mod_2 = npr % 2
    mpr_mod_2 = mpr % 2
    
    if n_mod_2 == m_mod_2 and npr_mod_2 == mpr_mod_2:
        k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q1Q2Q3Q4) 	

    elif n_mod_2 != m_mod_2:
        if npr_mod_2 == mpr_mod_2:
            k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q1Q4) + np.einsum('ij,ij', -u_xy_nodes, w_ij_Q2Q3)
        else:
            k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q2Q4) + np.einsum('ij,ij', -u_xy_nodes, w_ij_Q1Q3) 
            
    elif npr_mod_2 != mpr_mod_2:
        if n_mod_2 == m_mod_2:
            k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q3Q4) + np.einsum('ij,ij', -u_xy_nodes,  w_ij_Q1Q2)
        else:
            k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q2Q4) + np.einsum('ij,ij', -u_xy_nodes, w_ij_Q1Q3)
    
    return k_ROQ
    


def knmHG(couplings, q1, q2, surface_map=None, q1y=None, q2y=None, method="riemann"):
    if q1y == None:
        q1y = q1
        
    if q2y == None:
        q2y = q2
        
    assert q1.wavelength == q2.wavelength and q1y.wavelength == q2y.wavelength and q1y.wavelength == q1.wavelength
    
    couplings = np.array(couplings)
    
    a = couplings.size / 4.0
    
    if int(a) - a != 0:
        raise BasePyKatException("Iterator should be product of 4, each element of coupling array should be [n,m,n',m']")
        
    Axy = surface_map.z_xy(q1.wavelength)
    
    x = surface_map.x
    y = surface_map.y
    
    K = np.zeros((couplings.size/4,), dtype=np.complex128)
    
    it = np.nditer(couplings, flags=['refs_ok','f_index'])
    
    i = 0
    
    if method == "romhom":
        if surface_map == None:
            raise BasePyKatException("Using 'romhom' method requires a surface map to be specified")
            
        weights = surface_map.ROMWeights
        
        if weights == None:
            raise BasePyKatException("The ROM weights need to be generated for this map before use.")
        
        romcache = None#__gen_ROM_HG_knm_cache(weights, couplings, q1=q1, q2=q2, q1y=q1y, q2y=q2y)
    else:
        romcache = None
        weights = None
        
    while not it.finished:
        try:
            
            mode_in = [int(it.next()), int(it.next())]
            mode_out = [int(it.next()), int(it.next())]
            
            if method == "riemann":
                K[i] = riemann_HG_knm(x, y, mode_in, mode_out, q1=q1, q2=q2, q1y=q1y, q2y=q2y, Axy=Axy)
            elif method == "romhom":
                K[i] = ROM_HG_knm(weights, mode_in, mode_out, q1=q1, q2=q2, q1y=q1y, q2y=q2y, cache=romcache)
            else:
                raise BasePyKatException("method value '%s' not accepted" % method)
        
            i +=1
                 
        except StopIteration:
            break

    return K.reshape(couplings.shape[:-1])