import maps
import os.path
import numpy as np
import pykat
import collections
from itertools import combinations_with_replacement as combinations
from pykat.utilities.optics.gaussian_beams import beam_param, HG_beam
from scipy.linalg import inv
from math import factorial
from hermite import *
import math

EmpiricalInterpolant = collections.namedtuple('EmpiricalInterpolant', 'B nodes node_indices limits x')
ReducedBasis = collections.namedtuple('ReducedBasis', 'RB limits x')
ROMLimits = collections.namedtuple('ROMLimits', 'zmin zmax w0min w0max max_order')

class ROMWeights:
    
    def __init__(self, w_ij_Q1, w_ij_Q2, w_ij_Q3, w_ij_Q4, limits):
        self.w_ij_Q1 = w_ij_Q1
        self.w_ij_Q2 = w_ij_Q2
        self.w_ij_Q3 = w_ij_Q3
        self.w_ij_Q4 = w_ij_Q4
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

    return A_n1 * hermite(n1, np.sqrt(2.)*x / wz1) * np.exp(-1j*(2*math.pi/(1064e-9))* x**2 /(2.*q_z1))


def u_star_u(re_q1, re_q2, w0_1, w0_2, n1, n2, x):
    return u(re_q1, w0_1, n1, x) * u(re_q2, w0_2, n2, x).conjugate()

    
def makeReducedBasis(x, offset=np.array([0,0]), isModeMatched=True, tolerance = 1e-12, sigma = 1):
    
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

        q1 = beam_param(w0=params[i][0], z=params[i][2])

        if isModeMatched:
            q2 = q1
            n = int(params[i][2])
            m = int(params[i][3])
            w0_1 = params[i][0]
            w0_2 = w0_1
            re_q1 = params[i][1]
            re_q2 = re_q1
        else:
            q2 = beam_param(w0=params[i][1], z=params[i][3])            
            n = int(params[i][4])
            m = int(params[i][5])
            
        TS[IDx] = u_star_u(re_q1, re_q2, w0_1, w0_2, n, m, x)
        
        # normalize
        norm = np.sqrt(abs(np.vdot(dx * TS[IDx], TS[IDx])))
    
        if norm != 0:
            TS[IDx] /= norm 
            IDx += 1
        else:
            np.delete(TS, IDx)	

    #### Begin greedy: see Field et al. arXiv:1308.3565v2 #### 

    tolerance = 1e-12 # set maximum RB projection error

    sigma = 1 # projection error at 0th iteration

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
    
def makeEmpiricalInterpolant(RB):

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
    B = B_matrix(invV, e_x)

    return EmpiricalInterpolant(B=np.array(B), nodes=np.array(x_nodes), node_indices=np.array(node_indices),  limits=RB.limits, x=RB.x)
    
    
def makeWeights(smap, EI):
    
    # get full A_xy
    A_xy = smap.data
    
    hx = len(smap.x)/2
    hy = len(smap.y)/2
    fx = hx*2
    fy = hy*2
    
    # get A_xy in the four quadrants of the x-y plane
    A_xy_Q1 = A_xy[0:hx][0:hy, hy:fy]
    A_xy_Q2 = A_xy[hx:fx][0:hy, hy:fy]
    A_xy_Q3 =  A_xy[hx:fx][0:hy, 0:hy]
    A_xy_Q4 = A_xy[0:hx][0:hy, 0:hy]

    full_x = smap.x
    full_y = smap.y

    dx = full_x[1] - full_x[0]
    dy = full_y[1] - full_y[0]
   
    B = EI.B
    nodes_nv = EI.nodes 
    nodes_pv = - EI.nodes

    # make integration weights
    w_ij_Q1 = np.zeros((len(B), len(B)), dtype = complex)
    w_ij_Q2 = np.zeros((len(B), len(B)), dtype = complex)
    w_ij_Q3 = np.zeros((len(B), len(B)), dtype = complex)
    w_ij_Q4 = np.zeros((len(B), len(B)), dtype = complex)

    for i in range(len(B)):
        for j in range(len(B)):
            B_ij_Q1  = np.outer(B[i], B[j][::-1])		
            w_ij_Q1[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q1, A_xy_Q1)	

            B_ij_Q2 = np.outer(B[i][::-1], B[j][::-1])
            w_ij_Q2[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q2, A_xy_Q2)

            B_ij_Q3 = np.outer(B[i][::-1], B[j])
            w_ij_Q3[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q3, A_xy_Q3)

            B_ij_Q4 = np.outer(B[i], B[j])
            w_ij_Q4[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q4, A_xy_Q4)
            
    return ROMWeights(w_ij_Q1=w_ij_Q1, w_ij_Q2=w_ij_Q2, w_ij_Q3=w_ij_Q3, w_ij_Q4=w_ij_Q4, limits=EI.limits)
    
    
def ROMKnm(W, max_order, q1, q2, q1y=None, q2y=None):
    
    if q1y == None:
        q1y = q1
        
    if q2y == None:
        q2y = q2
        
    # get symmetric and anti-symmetric w_ij's 
    w_ij_Q1Q3 = W.w_ij_Q1 + W.w_ij_Q3
    w_ij_Q2Q4 = W.w_ij_Q2 + W.w_ij_Q4
    w_ij_Q1Q2 = W.w_ij_Q1 + W.w_ij_Q2
    w_ij_Q1Q4 = W.w_ij_Q1 + W.w_ij_Q4
    w_ij_Q2Q3 = W.w_ij_Q2 + W.w_ij_Q3
    w_ij_Q3Q4 = W.w_ij_Q3 + W.w_ij_Q4

    w_ij_Q1Q2Q3Q4 = W.w_ij_Q1 + W.w_ij_Q3 + W.w_ij_Q2 + W.w_ij_Q4

    num_fields = int((max_order + 1) * (max_order + 2) / 2);

    K_ROQ = np.array((num_fields, num_fields))

    for i in range(len(nm_pairs)):
        for j in range(len(nm_pairs)):

            # full quadrature
            n = nm_pairs[i][0]
            m = nm_pairs[i][1]
            npr = nm_pairs[j][0]
            mpr = nm_pairs[j][1]
            u_xy =  np.outer(u_star_u(re_q1, re_q2, w0_1, w0_2, n, m, full_x), u_star_u(re_q1, re_q2, w0_1, w0_2, npr, mpr, full_x))
            k = dx*dy*np.einsum('ij,ij', u_xy, A_xy)	
            
            # ROQ
            if nm_pairs[i][0] % 2 == 0 and nm_pairs[i][1] % 2 == 0 and nm_pairs[j][0] % 2 == 0 and nm_pairs[j][1] % 2 == 0:
                u_xy_nodes = np.outer(u_star_u(re_q1, re_q2, w0_1, w0_2, n, m, nodes_nv), u_star_u(re_q1, re_q2, w0_1, w0_2, npr, mpr, nodes_nv))
                k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q1Q2Q3Q4) 	
	
            elif nm_pairs[i][0] % 2 == 1 and nm_pairs[i][1] % 2 == 1 and nm_pairs[j][0] % 2 == 1 and nm_pairs[j][1] % 2 == 1:
                u_xy_nodes = np.outer(u_star_u(re_q1, re_q2, w0_1, w0_2, n, m, nodes_nv), u_star_u(re_q1, re_q2, w0_1, w0_2, npr, mpr, nodes_nv))
                k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q1Q2Q3Q4)

            elif nm_pairs[i][0] % 2 == 0 and nm_pairs[i][1] % 2 == 0 and nm_pairs[j][0] % 2 == 1 and nm_pairs[j][1] % 2 == 1:
                u_xy_nodes = np.outer(u_star_u(re_q1, re_q2, w0_1, w0_2, n, m, nodes_nv), u_star_u(re_q1, re_q2, w0_1, w0_2, npr, mpr, nodes_nv))
                k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q1Q2Q3Q4)

            elif nm_pairs[i][0] % 2 == 1 and nm_pairs[i][1] % 2 == 1 and nm_pairs[j][0] % 2 == 0 and nm_pairs[j][1] % 2 == 0:	
                u_xy_nodes = np.outer(u_star_u(re_q1, re_q2, w0_1, w0_2, n, m, nodes_nv), u_star_u(re_q1, re_q2, w0_1, w0_2, npr, mpr, nodes_nv))
                k_ROQ = np.einsum('ij,ij', u_xy_nodes, w_ij_Q1Q2Q3Q4)
                
            elif nm_pairs[i][0] % 2 == 0 and nm_pairs[i][1] % 2 == 1 or nm_pairs[i][0] % 2 == 1 and nm_pairs[i][1] % 2 == 0:
                u_x_nodes = u_star_u(re_q1, re_q2, w0_1, w0_2, n,m, nodes_nv)	
                u_y_nodes = u_star_u(re_q1, re_q2, w0_1, w0_2, npr, mpr, nodes_nv)   

                if nm_pairs[j][0] % 2 == 0 and nm_pairs[j][1] % 2 == 0 or nm_pairs[j][0] % 2 == 1 and nm_pairs[j][1] % 2 == 1:
                    u_xy_nodes_Q1Q4 = np.outer(u_x_nodes, u_y_nodes)
                    u_xy_nodes_Q2Q3 = - u_xy_nodes_Q1Q4
                    k_ROQ = np.einsum('ij,ij', u_xy_nodes_Q1Q4, w_ij_Q1Q4) + np.einsum('ij,ij', u_xy_nodes_Q2Q3, w_ij_Q2Q3)
                else:
                    u_xy_nodes_Q2Q4 = np.outer(u_x_nodes, u_y_nodes)
                    u_xy_nodes_Q1Q3 = - u_xy_nodes_Q2Q4
                    k_ROQ = np.einsum('ij,ij', u_xy_nodes_Q2Q4, w_ij_Q2Q4) + np.einsum('ij,ij', u_xy_nodes_Q1Q3, w_ij_Q1Q3) 

            elif nm_pairs[j][0] % 2 == 0 and nm_pairs[j][1] % 2 == 1 or nm_pairs[j][0] % 2 == 1 and nm_pairs[j][1] % 2 == 0:
                u_x_nodes = u_star_u(re_q1, re_q2, w0_1, w0_2, n, m, nodes_nv)   
                u_y_nodes = u_star_u(re_q1, re_q2, w0_1, w0_2, npr, mpr, nodes_nv)  
                
                if nm_pairs[i][0] % 2 == 0 and nm_pairs[i][1] % 2 == 0 or nm_pairs[i][0] % 2 == 1 and nm_pairs[i][1] % 2 == 1:
                    u_xy_nodes_Q3Q4 = np.outer(u_x_nodes, u_y_nodes)
                    u_xy_nodes_Q1Q2 = - u_xy_nodes_Q3Q4
                    k_ROQ = np.einsum('ij,ij', u_xy_nodes_Q3Q4, w_ij_Q3Q4) + np.einsum('ij,ij', u_xy_nodes_Q1Q2,  w_ij_Q1Q2)
                else:
                    u_xy_nodes_Q2Q4 = np.outer(u_x_nodes, u_y_nodes)
                    u_xy_nodes_Q1Q3 = - u_xy_nodes_Q2Q4
                    k_ROQ = np.einsum('ij,ij', u_xy_nodes_Q2Q4, w_ij_Q2Q4) + np.einsum('ij,ij', u_xy_nodes_Q1Q3, w_ij_Q1Q3)

