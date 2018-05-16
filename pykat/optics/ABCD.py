import numpy as np
from pykat.optics.gaussian_beams import BeamParam

def apply(ABCD, q1, n1, n2):
    """
    Applied the matrix ABCD to a complex beam parameter
    
    ABCD: Component matrix
    q1: input beam parameter
    n1: input refractive index
    n2: output refractive index
    """
    # make sure we have a complex here. We could have a beamParam object or a an actual complex value supplied
    #_q1 = complex(q1)
        
    return BeamParam(nr=n2, q=n2 * (ABCD[0,0] * q1/float(n1) + ABCD[0,1]) / (ABCD[1,0] * q1/float(n1) + ABCD[1,1]), wavelength=q1.wavelength)

def mirror_trans(n1, n2, Rc):
    return np.matrix([[1.0,0.0],[(n2-n1)/float(Rc),1.0]])

def mirror_refl(n, Rc):
    return np.matrix([[1.0,0.0],[-2.0*n/float(Rc),1.0]])

def bs_trans_x(n1, n2, Rcx, alpha, isDeg=True):
    if isDeg: alpha = np.deg2rad(alpha)
    
    alpha2=np.arcsin(float(n1)/float(n2)*np.sin(alpha))
    c1=np.cos(alpha)
    c2=np.cos(alpha2)
    delta_n_t=(n2*c2-n1*c1)/(c1*c2)
    return np.matrix([[c2/c1, 0.0], [delta_n_t/float(Rcx), c1/c2]])

def bs_trans_y(n1, n2, Rcy, alpha, isDeg=True):
    if isDeg: alpha = np.deg2rad(alpha)
    
    alpha2=np.arcsin(float(n1)/float(n2)*np.sin(alpha))
    c1=np.cos(alpha)
    c2=np.cos(alpha2)
    delta_n_s=n2*c2-n1*c1
    return np.matrix([[1.0, 0.0],[delta_n_s/float(Rcy),1.0]])
    
def bs_trans(n1, n2, Rcx, Rcy, alpha, isDeg=True):
    return (bs_trans_x(n1, n2, Rcx, alpha, isDeg), bs_trans_y(n1, n2, Rcy, alpha, isDeg))

def bs_refl_x(n, Rcx, alpha, isDeg=True):
    if isDeg: alpha = np.deg2rad(alpha)
    return np.matrix([[1.0, 0.0],[-2.0*n/(float(Rcx)*np.cos(alpha)),1.0]])
    
def bs_refl_y(n, Rcy, alpha, isDeg=True):
    if isDeg: alpha = np.deg2rad(alpha)
    return np.matrix([[1.0, 0.0],[-2.0*n*np.cos(alpha)/float(Rcy),1.0]])
        
def bs_refl(n1, Rcx, Rcy, alpha, isDeg=True):
    return (bs_refl_x(n1, Rcx, alpha, isDeg), bs_refl_y(n1, Rcy, alpha, isDeg))
    
def space(n,L):
    return np.matrix([[1.0, L/float(n)],[0.0,1.0]])
    
def lens(f):
    return np.matrix([[1.0, 0.0],[-1.0/f,1.0]])
