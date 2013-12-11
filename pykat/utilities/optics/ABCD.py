import numpy as np

def apply(ABCD, q1, n1, n2):
    return n2 * (ABCD[0,0] * q1/n1 + ABCD[0,1]) / (ABCD[1,0] * q1/n1 + ABCD[1,1])

def mirror_trans(n1, n2, Rc):
    return np.matrix([[1,0],[(n2-n1)/Rc,1]])
    
def mirror_refl(n, Rc):
    return np.matrix([[1,0],[-2*n/Rc,1]])
    
def space(n,L):
    return np.matrix([[1, L/n],[0,1]])
    
def lens(f):
    return np.matrix([[1, 0],[-1/f,1]])