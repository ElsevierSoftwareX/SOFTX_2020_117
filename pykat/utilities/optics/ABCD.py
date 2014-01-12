import numpy as np
from pykat.utilities.optics.gaussian_beams import gauss_param

def apply(ABCD, q1, n1, n2):
    return gauss_param(nr=n2, q=n2 * (ABCD[0,0] * q1/n1 + ABCD[0,1]) / (ABCD[1,0] * q1/n1 + ABCD[1,1]))

def mirror_trans(n1, n2, Rc):
    return np.matrix([[1,0],[(n2-n1)/Rc,1]])

def mirror_refl(n, Rc):
    return np.matrix([[1,0],[-2*n/Rc,1]])

def bs_trans(n1, n2, Rcx, Rcy, alpha)
    alpha2=math.asin(n1/n2*math.sin(alpha))
    c1=math.cos(alpha)
    c2=math.cos(alpha2)
    delta_n_t=(n2*c2-n1*c1)/(c1*c2)
    delta_n_s=n2*c2-n1*c1
    Mt= np.matrix([[c2/c1, 0],[delta_n_t/Rcx,c1/c2]])
    Ms= np.matrix([[1, 0],[delta_n_s/Rcx,1]])
    return(Mt, Ms)

def bs_refl(n1, n2, Rcx, Rcy, alpha)
    Mt= np.matrix([[1, 0],[-2*n1/(Rcx*math.cos(alpha)),1]])
    Ms= np.matrix([[1, 0],[-2*n1*math.cos(alpha)/Rcy,1]])
    return(Mt, Ms)
    
def space(n,L):
    return np.matrix([[1, L/n],[0,1]])
    
def lens(f):
    return np.matrix([[1, 0],[-1/f,1]])
