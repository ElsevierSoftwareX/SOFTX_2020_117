import math
import os.path
import pykat
import collections
    
from progressbar import ProgressBar, ETA, Percentage, Bar
from itertools import combinations_with_replacement as combinations
from pykat.utilities.optics.gaussian_beams import beam_param, HG_beam
from scipy.linalg import inv
from math import factorial
from hermite import *

import numpy as np

EmpiricalInterpolant = collections.namedtuple('EmpiricalInterpolant', 'B nodes node_indices limits x')
ReducedBasis = collections.namedtuple('ReducedBasis', 'RB limits x')
ROMLimits = collections.namedtuple('ROMLimits', 'zmin zmax w0min w0max max_order')

class ROMWeights:
    
    def __init__(self, w_ij_Q1, w_ij_Q2, w_ij_Q3, w_ij_Q4, EI, limits):
        self.w_ij_Q1 = w_ij_Q1
        self.w_ij_Q2 = w_ij_Q2
        self.w_ij_Q3 = w_ij_Q3
        self.w_ij_Q4 = w_ij_Q4
        
        self.EI = EI
        self.limits = limits
        
    def writeToFile(self, filename):
        f = open(filename, 'w+')
        
        f.write(repr(self.limits) + "\n")
        
        f.write(repr(self.w_ij_Q1.shape) + "\n")
        
        for v in self.w_ij_Q1.flatten():
            f.write("%s " % repr(complex(v))[1:-1])
        
        f.write("\n")
        
        f.write(repr(self.w_ij_Q2.shape) + "\n")
        
        for v in self.w_ij_Q2.flatten():
            f.write("%s " % repr(complex(v))[1:-1])
        
        f.write("\n")
        
        f.write(repr(self.w_ij_Q3.shape) + "\n")
        
        for v in self.w_ij_Q3.flatten():
            f.write("%s " % repr(complex(v))[1:-1])
        
        f.write("\n")
        
        f.write(repr(self.w_ij_Q4.shape) + "\n")
        
        for v in self.w_ij_Q4.flatten():
            f.write("%s " % repr(complex(v))[1:-1])
        
        f.write("\n")
            
        f.close()

def combs(a, r):
    """
    Return successive r-length combinations of elements in the array a.
    Should produce the same output as array(list(combinations(a, r))), but 
    faster.
    """
    a = np.asarray(a)
    dt = np.dtype([('', a.dtype)]*r)
    b = np.fromiter(combinations(a, r), dt)
    return b.view(a.dtype).reshape(-1, r)

def project_onto_basis(integration_weights, e, h, projections, proj_coefficients, idx):

    for j in range(len(h)):
        proj_coefficients[idx][j] = np.vdot(integration_weights* e[idx], h[j])
        projections[j] += proj_coefficients[idx][j]*e[idx]

    return projections
    
def B_matrix(invV, e):
    return np.inner(invV.T, e[0:(invV.shape[0])].T)

def emp_interp(B_matrix, func, indices):
    # B : RB matrix
    assert B_matrix.shape[0] == len(indices)
    interpolant = np.inner(func[indices].T, B_matrix.T)

    return interpolant 
   
    
def w(w0, im_q, re_q):
    return w0 * np.sqrt( 1 + (re_q / im_q)**2. )


def u(re_q1, w0_1, n1, x):
    
    im_q1 = np.pi*w0_1**2 / 1064e-9
    q_z1 = re_q1 + 1j*im_q1

    A_n1 = (2./np.pi)**(1./4.) * (1./((2.**n1)*factorial(n1)*w0_1))**(1./2.) * (im_q1 / q_z1)**(1./2.) * ( im_q1*q_z1.conjugate() / (-im_q1*q_z1)  )**(n1/2.) 

    wz1 = w(w0_1, im_q1, re_q1)

    return A_n1 * hermite(n1, np.sqrt(2.)*x / wz1) * np.exp(np.array(-1j*(2*math.pi/(1064e-9))* x**2 /(2.*q_z1)))


def u_star_u(re_q1, re_q2, w0_1, w0_2, n1, n2, x, x2=None):
    if x2 == None:
        x2 = x
        
    return u(re_q1, w0_1, n1, x) * u(re_q2, w0_2, n2, x2).conjugate()

    
def makeReducedBasis(x, isModeMatched=True, tolerance = 1e-12, sigma = 1, greedyfile=None):
    
    if greedyfile != None:
        greedypts = str(greedyfile)
    else:
        if isModeMatched:
            greedypts = 'matched20.txt'
        else:
            greedypts = 'mismatched20.txt'
    
    greedyfile = os.path.join(pykat.__path__[0],'utilities','greedypoints',greedypts)
    
    limits = np.loadtxt(greedyfile, usecols=(1,))[:5]
    
    romlimits = ROMLimits(w0min=limits[0], w0max=limits[1], zmin=limits[2], zmax=limits[3], max_order=limits[4])
    
    w_min = limits[0]
    w_max = limits[1]

    re_q_min = limits[2]
    re_q_max = limits[3]
    
    max_order = int(limits[4])
    
    indices = range(0, max_order+1)
    nm_pairs = combs(indices, 2)
    
    dx = x[1] - x[0]
    
    params = np.loadtxt(greedyfile, skiprows=5)

    TS_size = len(params) # training space of TS_size number of waveforms
    
    #### allocate memory for training space ####

    TS = np.zeros(TS_size*len(x), dtype = complex).reshape(TS_size, len(x)) # store training space in TS_size X len(x) array

    IDx = 0 #TS index 

    for i in range(len(params)):
        if isModeMatched:
            q1 = beam_param(w0=float(params[i][0]), z=float(params[i][1]))
            q2 = q1
            n = int(params[i][2])
            m = int(params[i][3])
        else:
            q1 = beam_param(w0=float(params[i][0]), z=float(params[i][2]))
            q2 = beam_param(w0=float(params[i][1]), z=float(params[i][3]))            
            n = int(params[i][4])
            m = int(params[i][5])
            
        w0_1 = q1.w0
        w0_2 = q2.w0
        re_q1 = q1.z
        re_q2 = q2.z
            
        TS[IDx] = u_star_u(re_q1, re_q2, w0_1, w0_2, n, m, x)
        
        # normalize
        norm = np.sqrt(abs(np.vdot(dx * TS[IDx], TS[IDx])))
    
        if norm != 0:
            TS[IDx] /= norm 
            IDx += 1
        else:
            np.delete(TS, IDx)	

    #### Begin greedy: see Field et al. arXiv:1308.3565v2 #### 

    RB_matrix = [TS[0]] # seed greedy algorithm (arbitrary)

    proj_coefficients = np.zeros(TS_size*TS_size, dtype = complex).reshape(TS_size, TS_size)
    projections = np.zeros(TS_size*len(x), dtype = complex).reshape(TS_size, len(x))

    iter = 0
    
    while sigma > tolerance:
    #for k in range(TS_size-1):
    	# go through training set and find worst projection error
    	projections = project_onto_basis(dx, RB_matrix, TS, projections, proj_coefficients, iter)
    	residual = TS - projections
        
    	projection_errors = [np.vdot(dx* residual[i], residual[i]) for i in range(len(residual))]
    	sigma = abs(max(projection_errors))
    	index = np.argmax(projection_errors) 
        
    	#Gram-Schmidt to get the next basis and normalize
    	next_basis = TS[index] - projections[index]
    	next_basis /= np.sqrt(abs(np.vdot(dx* next_basis, next_basis)))

    	RB_matrix.append(next_basis)		
        
    	iter += 1
        
        
    return ReducedBasis(RB=np.array(RB_matrix), limits=romlimits, x=x)
    
def makeEmpiricalInterpolant(RB, sort=False):

    e_x = RB.RB
    x = RB.x
    
    node_indices = []
    x_nodes = []
    V = np.zeros((len(e_x), len(e_x)), dtype = complex)
    
    # seed EIM algorithm
    node_indices.append( int(np.argmax( np.abs(e_x[0]) ))) 
    x_nodes.append(x[node_indices]) 
    
    for i in range(1, len(e_x)): 
        #build empirical interpolant for e_iter
        for j in range(len(node_indices)):  
            for k in range(len(node_indices)):  
                V[k][j] = e_x[j][node_indices[k]]  
            
        invV = inv(V[0:len(node_indices), 0:len(node_indices)])  
        B = B_matrix(invV, e_x) 
        
        interpolant = emp_interp(B, e_x[i], node_indices) 
        res = interpolant - e_x[i] 

        index = int(np.argmax(np.abs(res))) 
        node_indices.append(index) 
        x_nodes.append( x[index] )
    
    # make B matrix with all the indices
    for j in range(len(node_indices)):
        for k in range(len(node_indices)):
            V[k][j] = e_x[j][node_indices[k]]

    invV = inv(V[0:len(node_indices), 0:len(node_indices)])
    B = np.array(B_matrix(invV, e_x))
    
    node_indices = np.array(node_indices, dtype=np.int32)
    nodes = np.array(x_nodes, dtype=np.float64)
    
    if sort:
        print sort
        ix = np.argsort(nodes)
        nodes = nodes[ix]
        node_indices = node_indices[ix]
        B = B[ix, :]
    
    return EmpiricalInterpolant(B=B, nodes=nodes, node_indices=node_indices, limits=RB.limits, x=RB.x)
    

def makeWeights(smap, EI, verbose=True):
    
    # get full A_xy
    A_xy = smap.z_xy()
    
    xm = smap.x[smap.x < 0]
    xp = smap.x[smap.x > 0]
    
    ym = smap.y[smap.y < 0]
    yp = smap.y[smap.y > 0]
    
    Q1xy = np.ix_(smap.x < 0, smap.y > 0)
    Q2xy = np.ix_(smap.x > 0, smap.y > 0)
    Q3xy = np.ix_(smap.x > 0, smap.y < 0)
    Q4xy = np.ix_(smap.x < 0, smap.y < 0)
    
    # get A_xy in the four quadrants of the x-y plane
    A_xy_Q1 = A_xy[Q1xy]
    A_xy_Q2 = A_xy[Q2xy]
    A_xy_Q3 = A_xy[Q3xy]
    A_xy_Q4 = A_xy[Q4xy]
    
    # Not sure if there is a neater way to select the x=0,y=0 cross elements
    A_xy_0  = np.hstack([A_xy[:,smap.x==0].flatten(), A_xy[smap.y==0,:].flatten()])
    
    full_x = smap.x
    full_y = smap.y

    dx = full_x[1] - full_x[0]
    dy = full_y[1] - full_y[0]

    if verbose:
        count  = 4*len(EI["xm"].B) * len(EI["ym"].B)
        p = ProgressBar(maxval=count, widgets=["Computing weights: ", Percentage(), Bar(), ETA()])
    
    n = 0
    
    # make integration weights
    Bx = EI["xm"].B
    By = EI["ym"].B[:,::-1]
    w_ij_Q1 = np.zeros((len(Bx),len(By)), dtype = complex)
    
    for i in range(len(Bx)):
        for j in range(len(By)):
            B_ij_Q1 = np.outer(Bx[i], By[j])
            w_ij_Q1[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q1, A_xy_Q1)	
            
            if verbose:
                p.update(n)
                n+=1

    Bx = EI["xm"].B[:,::-1]
    By = EI["ym"].B[:,::-1]
    w_ij_Q2 = np.zeros((len(Bx),len(By)), dtype = complex)
    
    for i in range(len(Bx)):
        for j in range(len(By)):
            B_ij_Q2 = np.outer(Bx[i], By[j])
            w_ij_Q2[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q2, A_xy_Q2)
            
            if verbose:
                p.update(n)
                n+=1

    Bx = EI["xm"].B[:,::-1]
    By = EI["ym"].B
    w_ij_Q3 = np.zeros((len(Bx),len(By)), dtype = complex)
    
    for i in range(len(Bx)):
        for j in range(len(By)):
            B_ij_Q3 = np.outer(Bx[i], By[j])
            w_ij_Q3[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q3, A_xy_Q3)

            if verbose:
                p.update(n)
                n+=1

    Bx = EI["xm"].B
    By = EI["ym"].B
    w_ij_Q4 = np.zeros((len(Bx),len(By)), dtype = complex)
    
    for i in range(len(Bx)):
        for j in range(len(By)):
            B_ij_Q4 = np.outer(Bx[i], By[j])
            w_ij_Q4[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q4, A_xy_Q4)

            if verbose:
                p.update(n)
                n+=1
    if verbose:
        p.finish()
        
    return ROMWeights(w_ij_Q1=w_ij_Q1, w_ij_Q2=w_ij_Q2, w_ij_Q3=w_ij_Q3, w_ij_Q4=w_ij_Q4, EI=EI, limits=EI["xm"].limits)
    
    
    
    