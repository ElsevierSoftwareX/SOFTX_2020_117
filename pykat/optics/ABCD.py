import numpy as np
from pykat.optics.gaussian_beams import BeamParam

def apply(ABCD, q1, n1, n2):
    return BeamParam(nr=n2, q=n2 * (ABCD[0,0] * q1/float(n1) + ABCD[0,1]) / (ABCD[1,0] * q1/float(n1) + ABCD[1,1]))

def mirror_trans(n1, n2, Rc):
    return np.matrix([[1.0,0.0],[(n2-n1)/float(Rc),1.0]])

def mirror_refl(n, Rc):
    return np.matrix([[1.0,0.0],[-2.0*n/float(Rc),1.0]])

def bs_trans(n1, n2, Rcx, Rcy, alpha):
    alpha2=math.asin(n1/float(n2)*math.sin(alpha))
    c1=math.cos(alpha)
    c2=math.cos(alpha2)
    delta_n_t=(n2*c2-n1*c1)/(c1*c2)
    delta_n_s=n2*c2-n1*c1
    Mt= np.matrix([[c2/c1, 0.0],[delta_n_t/float(Rcx),c1/c2]])
    Ms= np.matrix([[1.0, 0.0],[delta_n_s/float(Rcy),1.0]])
    return(Mt, Ms)

def bs_refl(n1, n2, Rcx, Rcy, alpha):
    Mt= np.matrix([[1.0, 0.0],[-2.0*n1/(float(Rcx)*math.cos(alpha)),1.0]])
    Ms= np.matrix([[1.0, 0.0],[-2.0*n1*math.cos(alpha)/float(Rcy),1.0]])
    return(Mt, Ms)
    
def space(n,L):
    return np.matrix([[1.0, L/float(n)],[0.0,1.0]])
    
def lens(f):
    return np.matrix([[1.0, 0.0],[-1.0/f,1.0]])
