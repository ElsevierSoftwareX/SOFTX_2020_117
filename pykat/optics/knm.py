from itertools import combinations_with_replacement as combinations
from pykat.optics.gaussian_beams import BeamParam, HG_mode
from pykat.exceptions import BasePyKatException
from pykat.optics.romhom import u_star_u
from pykat.external.progressbar import ProgressBar, ETA, Percentage, Bar
from scipy.interpolate import interp2d
from scipy.integrate import dblquad
from pykat.optics.romhom import ROMWeights
from math import factorial
from pykat.math.hermite import hermite
from scipy.misc import comb
from scipy.integrate import newton_cotes
from pykat.math import newton_weights

import time
import pykat.optics.maps
import os.path
import numpy as np
import pykat
import collections
import math
import cmath

def makeCouplingMatrix(max_order, Neven=True, Nodd=True, Meven=True, Modd=True):
    max_order = int(max_order)
    c = []
    for n in range(0, max_order+1):
        for m in range(0, max_order+1):
            if n+m <= max_order:
                c.append([n,m])

    M = []

    for i in c:
        row = []
    
        for j in c:
            e = list(i)
            e.extend(j)
            
            if not Neven and (e[0]-e[2]) % 2 == 0: continue
            if not Nodd and (e[0]-e[2]) % 2 == 1: continue
            if not Meven and (e[1]-e[3]) % 2 == 0: continue
            if not Modd and (e[1]-e[3]) % 2 == 1: continue
            
            row.append(e)
        
        
        M.append(np.array(row).squeeze())
    
    return np.array(M)

def adaptive_knm(mode_in, mode_out, q1, q2, q1y=None, q2y=None, smap=None, delta=(0,0), params={}):
    
    if q1y is None:
        q1y = q1
        
    if q2y is None:
        q2y = q2
    
    if "epsabs" not in params: params["epsabs"] = 1e-6
    if "epsrel" not in params: params["epsrel"] = 1e-6
    if "usepolar" not in params: params["usepolar"] = False
        
    if len(mode_in) != 2 or len(mode_out) != 2:
        raise BasePyKatException("Both mode in and out should be a container with modes [n m]")
    
    Hg_in  = HG_mode(qx=q1, qy=q1y, n=mode_in[0], m=mode_in[1])
    Hg_out = HG_mode(qx=q2, qy=q2y, n=mode_out[0], m=mode_out[1])
    
    Nfuncs = []
    Nfuncs.append(0)
    
    if smap is not None:
        
        if not params["usepolar"]:
            xlims = (min(smap.x), max(smap.x))
            ylims = (min(smap.y), max(smap.y))
            
            def Rfunc(y,x):
                Nfuncs[-1] += len(x)
                return (Hg_in.Unm(x+delta[0], y+delta[1]) * smap.z_xy(x=x,y=y) * Hg_out.Unm(x, y).conjugate()).real
                
            def Ifunc(y,x):
                Nfuncs[-1] += len(x)
                return (Hg_in.Unm(x+delta[0], y+delta[1]) * smap.z_xy(x=x,y=y) * Hg_out.Unm(x, y).conjugate()).imag
            
        else:
            xlims = (0, 2*math.pi)
            ylims = (0, params["aperture"])
            
            def Rfunc(r, phi):
                Nfuncs[-1] += len(x)
                x = r*np.cos(phi)
                y = r*np.sin(phi)
                return (r * Hg_in.Unm(x, y) * smap.z_xy(x=x,y=y) * Hg_out.Unm(x, y).conjugate()).real
                
            def Ifunc(r, phi):
                Nfuncs[-1] += len(x)
                x = r*np.cos(phi)
                y = r*np.sin(phi)
                return (r * Hg_in.Unm(x, y) * smap.z_xy(x=x,y=y) * Hg_out.Unm(x, y).conjugate()).imag
            
    else:
        if not params["usepolar"]:
            _x = 4 * math.sqrt(1+max(mode_in[0],mode_in[1])) * q1.w
            _y = 4 * math.sqrt(1+max(mode_in[0],mode_in[1])) * q1y.w
        
            xlims = (-_x, _x)
            ylims = (-_y, _y)
        
            def Rfunc(y, x):
                Nfuncs[-1] += len(r)
                return (Hg_in.Unm(x+delta[0], y+delta[1]) * Hg_out.Unm(x, y).conjugate()).real
                
            def Ifunc(y,x):
                Nfuncs[-1] += len(r)
                return (Hg_in.Unm(x+delta[0], y+delta[1]) * Hg_out.Unm(x, y).conjugate()).imag
        else:
            xlims = (0, 2*math.pi)
            ylims = (0, params["aperture"])
            
            def Rfunc(r, phi):
                
                if hasattr(r, "__len__"):
                    Nfuncs[-1] += len(r)
                else:
                    Nfuncs[-1] += 1
                    
                x = r*np.cos(phi)
                y = r*np.sin(phi)
                return (r * Hg_in.Unm(x, y) * Hg_out.Unm(x, y).conjugate()).real
                
            def Ifunc(r, phi):
                if hasattr(r, "__len__"):
                    Nfuncs[-1] += len(r)
                else:
                    Nfuncs[-1] += 1
                    
                x = r*np.cos(phi)
                y = r*np.sin(phi)
                return (r * Hg_in.Unm(x, y) * Hg_out.Unm(x, y).conjugate()).imag
    
    R, errR = dblquad(Rfunc, xlims[0], xlims[1], lambda y: ylims[0], lambda y: ylims[1], epsabs=params["epsabs"], epsrel=params["epsrel"])
    I, errI = dblquad(Ifunc, xlims[0], xlims[1], lambda y: ylims[0], lambda y: ylims[1], epsabs=params["epsabs"], epsrel=params["epsrel"])
    
    params["Nfuncs"] = Nfuncs[0]
    params["errors"] = (errR, errI)
    
    return R + 1j * I
    
def riemann_HG_knm(x, y, mode_in, mode_out, q1, q2, q1y=None, q2y=None,
                     Axy=None, cache=None, delta=(0,0), params={}, newtonCotesOrder=0):

    if Axy is None:
        Axy == np.ones((len(x), len(y)))
    
    if q1y is None:
        q1y = q1
        
    if q2y is None:
        q2y = q2
        
    if len(mode_in) != 2 or len(mode_out) != 2:
        raise BasePyKatException("Both mode in and out should be a container with modes [n m]")        

    dx = abs(x[1] - x[0])
    dy = abs(y[1] - y[0])    
        
    if cache is None:
        Hg_in  = HG_mode(qx=q1, qy=q1y, n=mode_in[0], m=mode_in[1])
        Hg_out = HG_mode(qx=q2, qy=q2y, n=mode_out[0], m=mode_out[1])
        
        U1 = Hg_in.Unm(x+delta[0], y+delta[1])
        U2 = Hg_out.Unm(x,y).conjugate()

        if newtonCotesOrder > 0:
                
            W = newton_cotes(newtonCotesOrder, 1)[0]
            
            if newtonCotesOrder > 1:
                if (len(x) - len(W)) % newtonCotesOrder != 0:
                    raise ValueError("To use Newton-Cotes order {0} the number of data points in x must ensure: (N_x - ({0}+1)) mod {0} == 0".format(newtonCotesOrder) )

                if (len(y) - len(W)) % newtonCotesOrder != 0:
                    raise ValueError("To use Newton-Cotes order {0} the number of data points in y must ensure: (N_y - ({0}+1)) mod {0} == 0".format(newtonCotesOrder) )
                
            wx = np.zeros(x.shape, dtype=np.float64)    
            wy = np.zeros(y.shape, dtype=np.float64)
    
            N = len(W)

            for i in range(0, (len(wx)-1)/newtonCotesOrder): wx[(i*(N-1)):(i*(N-1)+N)] += W
            for i in range(0, (len(wy)-1)/newtonCotesOrder): wy[(i*(N-1)):(i*(N-1)+N)] += W
            
            Wxy = np.outer(wx, wy)
            
        if newtonCotesOrder == 0:
            return dx * dy * np.einsum('ij,ij', Axy, U1*U2)
        else:
            return dx * dy * np.einsum('ij,ij', Axy, U1*U2*Wxy)
    else:
        
        strx = "u1[%i,%i]" % (mode_in[0], mode_out[0])
        stry = "u2[%i,%i]" % (mode_in[1], mode_out[1])
        
        return dx * dy * np.einsum('ij,ij', Axy, np.outer(cache[strx], cache[stry]))



    
def __gen_riemann_knm_cache(x, y, couplings, q1, q2, q1y=None, q2y=None, delta=(0,0), params={}):
    if q1y is None:
        q1y = q1
        
    if q2y is None:
        q2y = q2
        
    #it = np.nditer(couplings, flags=['refs_ok','f_index'])
    
    cache = {}
    
    #while not it.finished:
    #    try:
    #        mode_in = [int(it.next()), int(it.next())]
    #        mode_out = [int(it.next()), int(it.next())]
    
    couplings = couplings.copy()
    couplings.resize(int(couplings.size/4), 4)
    
    for _ in couplings:
        mode_in = _[:2]
        mode_out = _[2:]        
        strx = "u1[%i,%i]" % (mode_in[0], mode_out[0])
        stry = "u2[%i,%i]" % (mode_in[1], mode_out[1])
        
        #Hg_in  = HG_beam(qx=q1, qy=q1y, n=mode_in[0], m=mode_in[1])
        #Hg_out = HG_beam(qx=q2, qy=q2y, n=mode_out[0], m=mode_out[1])

        if strx not in cache:
            cache[strx] = u_star_u(q1.z,   q2.z,  q1.w0,  q2.w0, mode_in[0], mode_out[0], x, x+delta[0])    
            #Hg_in.Un(x) * Hg_out.Un(x).conjugate()   
        
        if stry not in cache:
            cache[stry] = u_star_u(q1y.z,   q2y.z,  q1y.w0,  q2y.w0, mode_in[1], mode_out[1], y, y+delta[1])    
            #Hg_in.Um(y) * Hg_out.Um(y).conjugate()
            
    #    except StopIteration:
    #        break
    
    return cache
    
    
    
def __gen_ROM_HG_knm_cache(weights, couplings, q1, q2, q1y=None, q2y=None):

    if q1y is None:
        q1y = q1
        
    if q2y is None:
        q2y = q2
        
    
    
    cache = {}
    
    cache["w_ij_Q1Q3"] = weights.w_ij_Q1 + weights.w_ij_Q3
    cache["w_ij_Q2Q4"] = weights.w_ij_Q2 + weights.w_ij_Q4
    cache["w_ij_Q1Q2"] = weights.w_ij_Q1 + weights.w_ij_Q2
    cache["w_ij_Q1Q4"] = weights.w_ij_Q1 + weights.w_ij_Q4
    cache["w_ij_Q2Q3"] = weights.w_ij_Q2 + weights.w_ij_Q3
    cache["w_ij_Q3Q4"] = weights.w_ij_Q3 + weights.w_ij_Q4
    cache["w_ij_Q1Q2Q3Q4"] = weights.w_ij_Q1 + weights.w_ij_Q3 + weights.w_ij_Q2 + weights.w_ij_Q4

    couplings = couplings.copy()
    couplings.resize(int(couplings.size/4), 4)
    
    for _ in couplings:
        mode_in = _[:2]
        mode_out = _[2:]

        strx = "x[%i,%i]" % (mode_in[0], mode_out[0])
        stry = "y[%i,%i]" % (mode_in[1], mode_out[1])

        if strx not in cache:
            cache[strx] = u_star_u(q1.z,   q2.z,  q1.w0,  q2.w0, mode_in[0], mode_out[0], weights.EI["xm"].nodes)

        if stry not in cache:
            cache[stry] = u_star_u(q1y.z, q2y.z, q1y.w0, q2y.w0, mode_in[1], mode_out[1], weights.EI["ym"].nodes)
    
    # it = np.nditer(couplings, flags=['refs_ok','f_index'])
    # while not it.finished:
    #     try:
    #         mode_in = [int(it.next()), int(it.next())]
    #         mode_out = [int(it.next()), int(it.next())]
    #
    #         strx = "x[%i,%i]" % (mode_in[0], mode_out[0])
    #         stry = "y[%i,%i]" % (mode_in[1], mode_out[1])
    #
    #         if strx not in cache:
    #             cache[strx] = u_star_u(q1.z,   q2.z,  q1.w0,  q2.w0, mode_in[0], mode_out[0], weights.EI["xm"].nodes)
    #
    #         if stry not in cache:
    #             cache[stry] = u_star_u(q1y.z, q2y.z, q1y.w0, q2y.w0, mode_in[1], mode_out[1], weights.EI["ym"].nodes)
    #
    #     except StopIteration:
    #         break
    
    return cache



def ROM_HG_knm(weights, mode_in, mode_out, q1, q2, q1y=None, q2y=None, cache=None):
    if q1y is None:
        q1y = q1

    if q2y is None:
        q2y = q2
    
    # x modes
    n = mode_in[0]
    m = mode_out[0]

    # y modes
    npr = mode_in[1]
    mpr = mode_out[1]
    
    if isinstance(weights, ROMWeights):
        if cache is None:
            u_x_nodes = u_star_u(q1.z,   q2.z,  q1.w0,  q2.w0, n,     m,   weights.EIx.nodes)
            u_y_nodes = u_star_u(q1y.z,   q2y.z,  q1y.w0,  q2y.w0, npr, mpr,   weights.EIy.nodes)
        
            w_ij_Q1Q3 = weights.w_ij_Q1 + weights.w_ij_Q3
            w_ij_Q2Q4 = weights.w_ij_Q2 + weights.w_ij_Q4
            w_ij_Q1Q2 = weights.w_ij_Q1 + weights.w_ij_Q2
            w_ij_Q1Q4 = weights.w_ij_Q1 + weights.w_ij_Q4
            w_ij_Q2Q3 = weights.w_ij_Q2 + weights.w_ij_Q3
            w_ij_Q3Q4 = weights.w_ij_Q3 + weights.w_ij_Q4
            w_ij_Q1Q2Q3Q4 = weights.w_ij_Q1 + weights.w_ij_Q2 + weights.w_ij_Q3 + weights.w_ij_Q4
        
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
                k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q1Q4) - np.einsum('ij,ij', u_xy_nodes, w_ij_Q2Q3)
            else:
                k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q2Q4) - np.einsum('ij,ij', u_xy_nodes, w_ij_Q1Q3)

        elif npr_mod_2 != mpr_mod_2:
            if n_mod_2 == m_mod_2:
                k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q3Q4) - np.einsum('ij,ij', u_xy_nodes,  w_ij_Q1Q2)
            else:
                k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q2Q4) - np.einsum('ij,ij', u_xy_nodes, w_ij_Q1Q3)
    
    else:
        if cache is None:
            u_x_nodes = u_star_u(q1.z,   q2.z,  q1.w0,  q2.w0, n,     m,   weights.EIx.nodes)
            u_y_nodes = u_star_u(q1y.z,   q2y.z,  q1y.w0,  q2y.w0, npr, mpr,   weights.EIy.nodes)
    
            w_ij = weights.w_ij
        else:
            strx = "x[%i,%i]" % (mode_in[0], mode_out[0])
            stry = "y[%i,%i]" % (mode_in[1], mode_out[1])

            u_x_nodes = cache[strx]
            u_y_nodes = cache[stry]
    
        u_xy_nodes = np.outer(u_x_nodes, u_y_nodes)

        k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij)
         
    return k_ROQ

__fac_cache = []

def fac(n):
    global __fac_cache
    if len(__fac_cache) == 0:
        return math.factorial(int(n))
    else:
        return __fac_cache[n]

def m_1_pow(n):
    if n % 2 == 0:
        return 1
    else:
        return -1


def __Ss(u, _u, F, _F, d=0):
    r = 0
    
    for s in range(0, min(u,_u)+1):
        r += m_1_pow(s) * _F ** (u-s) * _F ** (_u-s) / (fac(2*s+d)*fac(u-s)*fac(_u-s))
        
    return r


def __S(m, _m, X, _X, F, _F, d=0):
    if m % 2 == 1:
        lim1 = int((m-1)/2)
    else:
        lim1 = int(m/2 )

    if _m % 2 == 1:
        lim2 = int((_m-1)/2)
    else:
        lim2 = int(_m/2)
    
    r = 0
    
    for u in range(0, lim1+1):
        for _u in range(0, lim2+1):
            r += m_1_pow(u) * _X**(m-2*u) * X**(_m-2*_u) / ( fac(m-2*u)*fac(_m-2*_u) )   * __Ss(u, _u, F, _F, d=d)
    
    return r
           

def __bayerhelms_kn(n, _n, q1, q2, gamma=0.0):
    
    K0 = (q1.zr - q2.zr)/q2.zr
    K2 = (q1.z - q2.z)/q2.zr
    K = (K0 + 1j*K2)/2.0
    
    Ktilde = abs(K / (1+K))

    if gamma != 0:
        a  = q2.zr * math.sin(gamma) / (cmath.sqrt(1+K.conjugate()) * q2.w0)

        _X = - a * (q2.z/q2.zr - 1j)
        X  = - a * (q2.z/q2.zr + 1j*(1+2*K.conjugate()))
        Ex = cmath.exp(-_X*X / 2.0)
    else:
        _X = 0.0
        X  = 0.0
        Ex = 1.0
    
    _F  = K / (2.0 * (1.0+K0))
    F = K.conjugate() / 2.0 

    Sg = __S(n, _n, X, _X, F, _F)

    if n > 0 and _n > 0:
        Su = __S(n-1, _n-1, X, _X, F, _F, 1)
    else:
        Su = 0
    
    b = m_1_pow(_n) * cmath.sqrt(fac(n) * fac(_n) * (1.0 + K.real)**(n+0.5) * (1.0 + K.conjugate()) ** (-(n+_n+1)))
    
    return b * Ex * (Sg - Su)


def bayerhelms_HG_knm(mode_in, mode_out, q1, q2, q1y=None, q2y=None, gamma=(0,0)):
    if q1y is None:
        q1y = q1

    if q2y is None:
        q2y = q2

    # x modes
    n = mode_in[0]
    _n = mode_out[0]

    # y modes
    m = mode_in[1]
    _m = mode_out[1]

    return __bayerhelms_kn(n,_n, q1, q2, 2*gamma[0]) * __bayerhelms_kn(m, _m, q1y, q2y, 2*gamma[1])

def __sqaure_knm_int(n, _n, R):
    # This uses the H_n(x) * H_m(x) product identity to reduce the overlap into
    # a sum of factorial and an integral of a single Hermite with a gaussian function
    # thus making it easier to solve
    expR = math.exp(-(R**2))
    S = 0
    
    for j in range(0, min(n, _n)+1):
        _j1 = _n + n - 2*j - 1
        
        if _j1+1 == 0:
            # for the zeroth order we just have the gaussian integral to solve
            L = math.sqrt(math.pi) * math.erf(R)    
        elif (_j1+1) % 2 == 1:
            # if the Hermite is odd then the integral is always 0 as its anti-symmetric
            L = 0
        else:
            L = 2 * hermite(_j1, 0) - expR * (hermite(_j1, R) - hermite(_j1, -R))
        
        I = 2**j * factorial(j) * comb(n, j) * comb(_n, j)
        
        S += I * L
                
    return S


def square_aperture_HG_knm(mode_in, mode_out, q, R):
    """
    Computes the coupling coefficients for a square aperture.
    """
    # x modes
    n = mode_in[0]
    _n = mode_out[0]

    # y modes
    m = mode_in[1]
    _m = mode_out[1]
    
    hg1 = HG_mode(q, n=n, m=m)
    hg2 = HG_mode(q, n=_n, m=_m)
        
    kx = hg1.constant_x * hg2.constant_x.conjugate()
    ky = hg1.constant_y * hg2.constant_y.conjugate()
     
    f = q.w / math.sqrt(2)
    R = R / (q.w / math.sqrt(2))
    
    kx *= f
    kx *= __sqaure_knm_int(n, _n, R)
    
    ky *= f
    ky *= __sqaure_knm_int(m, _m, R)
    
    return kx * ky



def knmHG(couplings, q1, q2, surface_map=None, q1y=None, q2y=None, method="riemann", verbose=False, profile=False, gamma=(0,0), delta=(0,0), params={}):
    if q1y is None:
        q1y = q1
        
    if q2y is None:
        q2y = q2
        
    assert q1.wavelength == q2.wavelength and q1y.wavelength == q2y.wavelength and q1y.wavelength == q1.wavelength
    
    couplings = np.array(couplings)
    
    a = couplings.size / 4.0
    
    if int(a) - a != 0:
        raise BasePyKatException("Iterator should be product of 4, each element of coupling array should be [n,m,n',m']")
    
    maxtem = 0
    c = couplings.flatten()
    
    for i in range(0, int(c.size/2)):
        maxtem = max(sum(c[i*2:(i*2+2)]), maxtem)
    
    global __fac_cache
    
    for n in range(0, maxtem+1):
        __fac_cache.append(math.factorial(n))
    
    if surface_map is not None:  
        Axy = surface_map.z_xy(wavelength=q1.wavelength)
    
        x = surface_map.x
        y = surface_map.y
    
    K = np.zeros((couplings.size/4,), dtype=np.complex128)
    
    #it = np.nditer(couplings, flags=['refs_ok','f_index'])
    
    i = 0
    
    if profile:
        t0 = time.time()
        
    if method == "romhom":
        if surface_map is None:
            raise BasePyKatException("Using 'romhom' method requires a surface map to be specified")
            
        weights = surface_map.ROMWeights
        
        if weights is None:
            raise BasePyKatException("The ROM weights need to be generated for this map before use.")

        cache = __gen_ROM_HG_knm_cache(weights, couplings, q1=q1, q2=q2, q1y=q1y, q2y=q2y)
        
    elif method == "riemann":
        if surface_map is None:
            raise BasePyKatException("Using 'riemann' method requires a surface map to be specified")
            
        cache = __gen_riemann_knm_cache(x, y, couplings, q1, q2, q1y=None, q2y=None, delta=delta)
    else:
        cache = None
        weights = None
    
    if profile:
        cache_time = time.time() - t0
        Ktime = np.zeros((couplings.size/4,), dtype=np.float64)
    
    if verbose:
        p = ProgressBar(maxval=couplings.size, widgets=["Knm (%s): " % method, Percentage(), Bar(), ETA()])
    
    _couplings = couplings.copy()
    _couplings.resize(int(_couplings.size/4), 4)
    
    for _ in _couplings:
        mode_in = _[:2]
        mode_out = _[2:]

        if profile:
            t0 = time.time()
                
        
        
        if method == "riemann":
            K[i] = riemann_HG_knm(x, y, mode_in, mode_out, q1=q1, q2=q2, q1y=q1y, q2y=q2y, Axy=Axy, cache=cache, delta=delta)
        elif method == "romhom":
            K[i] = ROM_HG_knm(weights, mode_in, mode_out, q1=q1, q2=q2, q1y=q1y, q2y=q2y, cache=cache)
        elif method == "bayerhelms":
            K[i] = bayerhelms_HG_knm(mode_in, mode_out, q1=q1, q2=q2, q1y=q1y, q2y=q2y, gamma=gamma)
        elif method == "adaptive":
            K[i] = adaptive_knm(mode_in, mode_out, q1=q1, q2=q2, q1y=q1y, q2y=q2y, smap=surface_map, delta=delta, params=params)
        else:
            raise BasePyKatException("method value '%s' not accepted" % method)
        
        if profile:
            Ktime[i] = time.time() - t0
        
        i +=1
        
        if verbose:
            p.update(i*4)

    if profile:
        return K.reshape(couplings.shape[:-1]), Ktime.reshape(couplings.shape[:-1]), cache_time
    else:
        return K.reshape(couplings.shape[:-1])


def plot_knm_matrix(couplings, knm, cmap=None, show=True):
    import pylab as plt
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.pcolormesh(abs(knm), cmap=cmap)
    fig.colorbar(cax)
    
    numrows, numcols = knm.shape
    
    c = couplings[:, 0, :2]
    c_ = []

    for d in c:
        c_.append("[%i,%i]"%(d[0], d[1]))
    
    A = np.arange(1, len(c)+1)-0.5
    ax.set_xticks(A)
    ax.set_yticks(A)
    ax.set_xticklabels(c_)
    ax.set_yticklabels(c_)
    
    ax.set_xlim(None, max(A)+0.5)
    ax.set_ylim(None, max(A)+0.5)
    
    def format_coord(x, y):
        col = int(np.floor(x))
        row = int(np.floor(y))
        
        if col>=0 and col<numcols and row>=0 and row<numrows:
            z = knm[row,col]
            return 'x=%s, y=%s, z=%1.4f' % (c_[col], c_[row], z)
        
        return None

    ax.format_coord = format_coord

    fig.tight_layout()
    
    if show: plt.show()
    
    return fig
