import math
import os.path
import pykat
import collections
import numpy as np
import multiprocessing
import h5py
import time
import datetime
import pickle
import itertools

from copy import copy
from pykat.external.progressbar import ProgressBar, ETA, Percentage, Bar
from itertools import combinations_with_replacement as combinations
from pykat.optics.gaussian_beams import beam_param, HG_beam
from scipy.linalg import inv
from math import factorial
from pykat.maths.hermite import *
from pykat.maths import newton_weights
from scipy.integrate import newton_cotes
from multiprocessing import Process, Queue, Array, Value, Event


EmpiricalInterpolant = collections.namedtuple('EmpiricalInterpolant', 'B nodes node_indices limits x')
ReducedBasis = collections.namedtuple('ReducedBasis', 'RB limits x')
ROMLimits = collections.namedtuple('ROMLimits', 'zmin zmax w0min w0max R mapSamples newtonCotesOrder max_order')
                       
class ROMWeights:
    
    def __init__(self, w_ij_Q1, w_ij_Q2, w_ij_Q3, w_ij_Q4, EIx, EIy):
        self.w_ij_Q1 = w_ij_Q1
        self.w_ij_Q2 = w_ij_Q2
        self.w_ij_Q3 = w_ij_Q3
        self.w_ij_Q4 = w_ij_Q4
        
        self.EIx = EIx
        self.EIy = EIy
        
    def writeToFile(self, filename):
        """
        Writes this map's ROM weights to a file
        that can be used with Finesse. The filename
        is appended with '.rom' internally.
        """
        f = open(filename + ".rom", 'w+')
        
        f.write("zmin=%16.16e\n" % self.EIx.limits.zmin)
        f.write("zmax=%16.16e\n" % self.EIx.limits.zmax)
        f.write("w0min=%16.16e\n" % self.EIx.limits.w0min)
        f.write("w0max=%16.16e\n" % self.EIx.limits.w0max)
        f.write("maxorder=%i\n" % self.EIx.limits.max_order)
        
        f.write("xnodes=%i\n" % len(self.EIx.nodes))
        
        for v in self.EIx.nodes.flatten():
            f.write("%s\n" % repr(float(v)))
        
        f.write("ynodes=%i\n" % len(self.EIy.nodes))
        
        for v in self.EIy.nodes.flatten():
            f.write("%s\n" % repr(float(v)))
            
        f.write(repr(self.w_ij_Q1.shape) + "\n")
        
        for v in self.w_ij_Q1.flatten():
            f.write("%s\n" % repr(complex(v))[1:-1])
        
        f.write(repr(self.w_ij_Q2.shape) + "\n")
        
        for v in self.w_ij_Q2.flatten():
            f.write("%s\n" % repr(complex(v))[1:-1])
        
        f.write(repr(self.w_ij_Q3.shape) + "\n")
        
        for v in self.w_ij_Q3.flatten():
            f.write("%s\n" % repr(complex(v))[1:-1])
        
        f.write(repr(self.w_ij_Q4.shape) + "\n")
        
        for v in self.w_ij_Q4.flatten():
            f.write("%s\n" % repr(complex(v))[1:-1])
            
        f.close()

class ROMWeightsFull:
    
    def __init__(self, w_ij, EI, limits):
        self.w_ij = w_ij
        
        self.EI = EI
        self.limits = limits
        
    def writeToFile(self, filename):
        """
        Writes this map's ROM weights to a file
        that can be used with Finesse. The filename
        is appended with '.rom' internally.
        """
        f = open(filename + ".rom", 'w+')
        
        f.write("zmin=%16.16e\n" % self.limits.zmin)
        f.write("zmax=%16.16e\n" % self.limits.zmax)
        f.write("w0min=%16.16e\n" % self.limits.w0min)
        f.write("w0max=%16.16e\n" % self.limits.w0max)
        f.write("maxorder=%i\n" % self.limits.max_order)
        
        f.write("xnodes=%i\n" % len(self.EI["x"].nodes))
        
        for v in self.EI["x"].nodes.flatten():
            f.write("%s\n" % repr(float(v)))
        
        f.write("ynodes=%i\n" % len(self.EI["y"].nodes))
        
        for v in self.EI["y"].nodes.flatten():
            f.write("%s\n" % repr(float(v)))
            
        f.write(repr(self.w_ij.shape) + "\n")
        
        for v in self.w_ij.flatten():
            f.write("%s\n" % repr(complex(v))[1:-1])
        
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
    if x2 is None:
        x2 = x
        
    return u(re_q1, w0_1, n1, x) * u(re_q2, w0_2, n2, x2).conjugate()

def u_star_u_mm(z, w0, n1, n2, x):
    return u(z, w0, n1, x) * u(z, w0, n2, x).conjugate()
    
################################################################################################
################################################################################################
################################################################################################

def makeReducedBasis(x, isModeMatched=True, tolerance = 1e-12, sigma = 1, greedyfile=None):
    
    if greedyfile is not None:
        greedypts = str(greedyfile)
    else:
        if isModeMatched:
            greedypts = 'matched20.txt'
        else:
            greedypts = 'mismatched20.txt'
    
    greedyfile = os.path.join(pykat.__path__[0],'optics','greedypoints',greedypts)
    
    # Two potential different formats
    try:
        limits = np.loadtxt(greedyfile, usecols=(3,))[:5]
    except IndexError:
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
    

def makeWeights(smap, EI, verbose=True, useSymmetry=True, newtonCotesOrder=1):
    
    if useSymmetry:
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
            count  = 4*len(EI["x"].B) * len(EI["y"].B)
            p = ProgressBar(maxval=count, widgets=["Computing weights: ", Percentage(), Bar(), ETA()])
    
        n = 0
    
        # make integration weights
        Bx = EI["x"].B
        By = EI["y"].B[:,::-1]
        w_ij_Q1 = np.zeros((len(Bx),len(By)), dtype = complex)

        from pykat.optics.knm import newton_weights
        
        wx = newton_weights(xm, newtonCotesOrder)
        wy = newton_weights(ym, newtonCotesOrder)
        
        for i in range(len(Bx)):
            for j in range(len(By)):
                B_ij_Q1 = np.outer(wx*Bx[i], wy*By[j])
                w_ij_Q1[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q1, A_xy_Q1)	
            
                if verbose:
                    p.update(n)
                    n+=1

        Bx = EI["x"].B[:,::-1]
        By = EI["y"].B[:,::-1]
        w_ij_Q2 = np.zeros((len(Bx),len(By)), dtype = complex)
    
        for i in range(len(Bx)):
            for j in range(len(By)):
                B_ij_Q2 = np.outer(wx*Bx[i], wy*By[j])
                w_ij_Q2[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q2, A_xy_Q2)
            
                if verbose:
                    p.update(n)
                    n+=1

        Bx = EI["x"].B[:,::-1]
        By = EI["y"].B
        w_ij_Q3 = np.zeros((len(Bx),len(By)), dtype = complex)
    
        for i in range(len(Bx)):
            for j in range(len(By)):
                B_ij_Q3 = np.outer(wx*Bx[i], wy*By[j])
                w_ij_Q3[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q3, A_xy_Q3)

                if verbose:
                    p.update(n)
                    n+=1

        Bx = EI["x"].B
        By = EI["y"].B
        w_ij_Q4 = np.zeros((len(Bx),len(By)), dtype = complex)
    
        for i in range(len(Bx)):
            for j in range(len(By)):
                B_ij_Q4 = np.outer(wx*Bx[i], wy*By[j])
                w_ij_Q4[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q4, A_xy_Q4)

                if verbose:
                    p.update(n)
                    n+=1
        if verbose:
            p.finish()
        
        return ROMWeights(w_ij_Q1=w_ij_Q1, w_ij_Q2=w_ij_Q2, w_ij_Q3=w_ij_Q3, w_ij_Q4=w_ij_Q4, EI=EI, limits=EI["x"].limits)
        
    else:
        # get full A_xy
        A_xy = smap.z_xy()

        full_x = smap.x
        full_y = smap.y

        dx = full_x[1] - full_x[0]
        dy = full_y[1] - full_y[0]

        if verbose:
            count  = len(EI["x"].B) * len(EI["y"].B)
            p = ProgressBar(maxval=count, widgets=["Computing weights: ", Percentage(), Bar(), ETA()])

        n = 0

        # make integration weights
        Bx = EI["x"].B
        By = EI["y"].B
        
        w_ij = np.zeros((len(Bx), len(By)), dtype = complex)

        for i in range(len(Bx)):
            for j in range(len(By)):
                B_ij = np.outer(Bx[i], By[j])
                w_ij[i][j] = dx*dy*np.einsum('ij,ij', B_ij, A_xy)	
    
                if verbose:
                    p.update(n)
                    n+=1
                    
        if verbose:
            p.finish()

        return ROMWeightsFull(w_ij=w_ij, EI=EI, limits=EI["x"].limits)


###################################################################################################
###################################################################################################
###################################################################################################
# !!! New ROM code below that doesn't need supercomputer
###################################################################################################
###################################################################################################
###################################################################################################
def _compute_TS(queue, oqueue, x, w):
    while True:
        msg = queue.get()
        
        if msg is None:
            break
        else:
            tmp = u_star_u_mm(msg[0], msg[1], msg[2], msg[3], x)
            # includes normalisation with quadrature rule weights
            norm = np.sqrt(1/(abs(np.vdot(w*tmp,tmp))))
            oqueue.put((msg, tmp*norm))


def _write_TS(queue, filename, tssize):
    hfile = h5py.File("%s.h5" % filename, 'a') 
    
    i = 0
    
    try:
        while True:
            msg = queue.get()
            
            if msg is None:
                break
            else:
                # Dump each TS into a group
                key = 'TS/%i' % msg[0][4]
                
                hfile[key+"/data"] = msg[1]
                hfile[key+"/z"]  = msg[0][0]
                hfile[key+"/w0"] = msg[0][1]
                hfile[key+"/n1"] = msg[0][2]
                hfile[key+"/n2"] = msg[0][3]
                
                i += 1
                
                if i % 1000 == 0:
                    print i/float(tssize)
                    hfile.flush()
    finally:
        hfile.close()
        
def CreateTrainingSetHDF5(filename, maxOrder, z, w0, R, halfMapSamples, newtonCotesOrder=1, NProcesses=1):
    """
    newtonCotesOrder: Order of integration to use
                        0 - Midpoint
                        1 - Trapezoid
                        2 - Simpsons
                        Or higher orders
    """
    
    iq = Queue()
    oq = Queue()
        
    Ns = halfMapSamples
    
    h = R / float(Ns-1) # step size
    
    # Close newton-cotes quadrature goes right upto the boundary
    # unlike previous midpoint rule.
    x = np.linspace(-R, 0, Ns, dtype=np.float64)
    
    w = newton_weights(x, newtonCotesOrder)

    nModes = 0
    
    for n in xrange(0, maxOrder+1):
        for m in xrange(0, maxOrder+1):
            if n+m <= maxOrder and n <= m:
                nModes += 1
            
    tssize = len(w0) * len(z) * nModes
    
    hfile = h5py.File("%s.h5" % filename, 'w')
     
    hfile['TSSize'] = tssize
    hfile['x'] = x
    hfile['zRange'] = (min(z), max(z))
    hfile['w0Range'] = (min(w0), max(w0))
    hfile['R'] = R
    hfile['halfMapSamples'] = halfMapSamples
    hfile['maxOrder'] = maxOrder
    hfile['newtonCotesOrder'] = newtonCotesOrder
    hfile['weights'] = w
    
    hfile.close() # make sure it's closed before
    
    print "Starting processes..."
    # Have a bunch of processes doing the computation and one doing the writing
    iprocesses = [Process(target=_compute_TS, name="irom%i" % i, args=(iq, oq, x, w)) for i in range(NProcesses)]
    oprocess = Process(target=_write_TS, name="orom%i" % i, args=(oq, filename, tssize))
    
    oprocess.start()
    
    try:
        for P in iprocesses:
            P.start()
        
        print "Filling queue..."
        curr = 0
        for n in xrange(0, maxOrder+1):
            for m in xrange(0, maxOrder+1):
                if n+m <= maxOrder and n <= m:
                    for _z in z:
                        for _w0 in w0:
                            iq.put((_z, _w0, n, m, curr))
                            # Can't use Queue.qsize() on OSX so count normally...
                            curr += 1
        
        for P in iprocesses:
            iq.put(None) # put a None for each thread to catch to end
            
        for P in iprocesses:
            P.join()
    
        # Put none to stop output process
        oq.put(None)
        
        oprocess.join()
        
    except:
        print("Exception occurred")
        
        for P in iprocesses:
            P.terminate()
        
        oprocess.terminate()
    
    print("Completed training set creation")
    print("Data written to %s.h5" % filename)
    
    
def _worker_ROM(hdf5Filename, job_queue, result_err, result_idx, event):
    # h5py drivers: 'core', 'sec2', 'stdio', 'mpio'
    # Need to use something ot her than sec2, the default on OSX,
    # as it doesn't play nice with multiple processes reading files
    with h5py.File("%s.h5" % hdf5Filename, driver="core", mode='r') as file:
    
        TS = file["TS"]
        
        while True:
            
            msg = job_queue.get()
            
            if msg is None:
                break
            else:
                TSidx, B, EI_indices = msg
                TSidx = np.array(TSidx)
                
                max_err = []
                
                for ll in TSidx:
                    a = TS['%i/data' % ll].value
                    
                    _err = np.max(np.abs(a - emp_interp(B, a, EI_indices)))
                        
                    max_err.append(_err)
                    
                result_err.value = np.max(max_err)
                result_idx.value = TSidx[np.argmax(max_err)]
                
                event.set()

def MakeROMFromHDF5(hdf5Filename, greedyFilename=None, EIFilename=None, tol=1e-10, NProcesses=1, maxRBsize=50):
    start = time.time()
    
    #### Start reading TS file ####
    TSdata = h5py.File("%s.h5" % hdf5Filename, 'r') 
    TS = TSdata['TS']		
    
    quadratureWeights = TSdata['weights'][...]
    x = TSdata['x'][...]
    maxRBsize = maxRBsize
    TSsize = TSdata['TSSize'][...]
    
    #### Set up stuff for greedy #### 
    tol = tol
    rb_errors = []
    x_nodes = []
    
    # Initial RB to seed with
    next_RB_index = 0  
    #EI_indices = [next_RB_index]
    #x_nodes.append(x[EI_indices])
	
    RB0 = TS['%i/data' % next_RB_index][...]
    #RB0 /= RB0[EI_indices[0]]
    EI_indices = [np.argmax(RB0)]
    x_nodes.extend(x[EI_indices])
    RB0 /= RB0[EI_indices[0]]		
    RB_matrix = [RB0] 

    V = np.zeros(((maxRBsize), (maxRBsize)), dtype=complex)

    V[0][0] = RB_matrix[0][EI_indices[0]]
    invV = inv(V[0:len(EI_indices), 0:len(EI_indices)])
    B = B_matrix(invV, np.array(RB_matrix))	
    
    RBs = []
    
    if NProcesses > 1:
        queue = Queue()
        locks = [Event() for l in range(NProcesses)]
        
        result_err = [Value('d', np.inf) for l in range(NProcesses)]
        result_idx = [Value('i', -1)     for l in range(NProcesses)]
        
        Names = range(NProcesses)
        procs = [Process(name="process_%i" % l[0], target=_worker_ROM, 
                         args=(hdf5Filename, queue, l[1], l[2], l[3])) for l in itertools.izip(Names, result_err, result_idx, locks)]

        max_res = np.zeros((NProcesses), dtype='d')
        max_idx = np.zeros((NProcesses), dtype='i')
        
        for P in procs: P.start()
    
    dstr = datetime.datetime.strftime(datetime.datetime.now(), "%d%m%Y_%H%M%S")
    
    if greedyFilename is None:
        greedyFilename = "GreedyPoints_%s" % dstr
    
    greedyFilename += ".dat"
    
    limits = ROMLimits(zmin=min(TSdata['zRange'].value),
                       zmax=max(TSdata['zRange'].value),
                       w0min=min(TSdata['w0Range'].value),
                       w0max=max(TSdata['w0Range'].value),
                       R=TSdata['R'].value,
                       mapSamples=TSdata['halfMapSamples'].value,
                       newtonCotesOrder=int(TSdata['newtonCotesOrder'].value),
                       max_order=int(TSdata['maxOrder'].value))
                       
    with open(greedyFilename, "w") as f:
        f.write("min w0 = %15.15e\n" % limits.zmin)
        f.write("max w0 = %15.15e\n" % limits.zmax)
        f.write("min z  = %15.15e\n" % limits.w0min)
        f.write("min z  = %15.15e\n" % limits.w0max)
        f.write("R      = %15.15e\n" % limits.R)
        f.write("max order   = %i\n" % limits.max_order)
        f.write("NC order    = %i\n" % limits.newtonCotesOrder)
        f.write("half map samples = %i\n" % limits.mapSamples)
        
        # write initial RB
        _TS = TS[str(next_RB_index)]
        f.write("%15.15e %15.15e %i %i\n" % (_TS["z"].value, _TS["w0"].value, _TS["n1"].value, _TS["n2"].value))
    
        for k in range(1, maxRBsize): 
        
            if NProcesses == 1:
                max_res = []

                TSidx = np.array(range(TSsize))

                for ll in TSidx:
                    a = TS['%i/data' % ll].value
                    max_res.append(np.max(a - emp_interp(B, a, EI_indices)))

                worst_error = max(np.abs(max_res))
                next_RB_index = np.argmax(max_res)
            
            else:
            
                TSs = [range(TSsize)[i::NProcesses] for i in range(NProcesses)]
            
                for l in TSs: queue.put((l, B, EI_indices))
            
                end_locks = copy(locks)
            
                while len(end_locks) > 0:
                    for e in end_locks:
                        if e.wait():
                            end_locks.remove(e)
            
                for e in locks: e.clear()
                
                for i in range(NProcesses):
                    max_res[i] = result_err[i].value
                    max_idx[i] = result_idx[i].value
            
                worst_error = max(np.abs(max_res))
                next_RB_index = max_idx[np.argmax(max_res)]
        
            if worst_error <= tol:
                print "Final basis size = %d, Final error = %e, Tolerance=%e" % (k, worst_error, tol) 
                break

            print "worst error = %e at %i on iteration %d" % (worst_error, next_RB_index, k)			
        
            epsilon = TS['%i/data' % next_RB_index].value
            res   = epsilon - emp_interp(B, epsilon, EI_indices)
            index = np.argmax(res)
            EI_indices.append(index)
            x_nodes.append(TSdata['x'][index])
            RB_matrix.append(res/max(res))

            for l in range(len(EI_indices)): # Part of (5) of Algorithm 2: making V_{ij} 
                for m in range(len(EI_indices)): # Part of (5) of Algorithm 2: making V_{ij} 
                    V[m][l] = RB_matrix[l][EI_indices[m]] # Part of (5) of Algorithm 2: making V_{ij}

            invV = inv(V[0:len(EI_indices), 0:len(EI_indices)])
            B = B_matrix(invV, np.array(RB_matrix))
            
            _TS = TS[str(next_RB_index)]
            f.write("%15.15e %15.15e %i %i\n" % (_TS["w0"].value, _TS["z"].value, _TS["n1"].value, _TS["n2"].value))
                
    print time.time() - start, "Seconds"
    
    if NProcesses > 1:
        for P in procs: P.terminate()

    TSdata.close()

    greedyFilenameBase = os.path.splitext(greedyFilename)[0]
    
    print "Writing to %s" % greedyFilename
                       
    EI = EmpiricalInterpolant(B=np.matrix(B).real,
                              nodes=np.array(x_nodes).squeeze(),
                              node_indices=np.array(EI_indices).squeeze(),
                              limits=limits,
                              x=x.squeeze())
    
    if EIFilename is not None:
        with open("%s.p" % EIFilename, 'wb') as f:
            pickle.dump(EI, f)
        
        print "Writing to %s.p" % EIFilename
                              
    return EI
    
    
def makeWeightsNew(smap, EIxFilename, EIyFilename=None, verbose=True, newtonCotesOrder=1):
    with open("%s" % EIxFilename, 'rb') as f:
        EIx = pickle.load(f)
        
    if EIyFilename is None:
        EIy = EIx
    else:
        with open("%s" % EIyFilename, 'rb') as f:
            EIy = pickle.load(f)

    # get full A_xy
    A_xy = smap.z_xy()

    xm = smap.x[smap.x <= 0]
    xp = smap.x[smap.x >= 0]
    ym = smap.y[smap.y <= 0]
    yp = smap.y[smap.y >= 0]

    Q1xy = np.ix_(smap.x <= 0, smap.y >= 0)
    Q2xy = np.ix_(smap.x >= 0, smap.y >= 0)
    Q3xy = np.ix_(smap.x >= 0, smap.y <= 0)
    Q4xy = np.ix_(smap.x <= 0, smap.y <= 0)

    # get A_xy in the four quadrants of the x-y plane
    A_xy_Q1 = A_xy[Q1xy]
    A_xy_Q2 = A_xy[Q2xy]
    A_xy_Q3 = A_xy[Q3xy]
    A_xy_Q4 = A_xy[Q4xy]

    full_x = smap.x
    full_y = smap.y
    
    dx = full_x[1] - full_x[0]
    dy = full_y[1] - full_y[0]

    if verbose:
        count  = 4*len(EIx.B) * len(EIy.B)
        p = ProgressBar(maxval=count, widgets=["Computing weights: ", Percentage(), Bar(), ETA()])

    n = 0

    # make integration weights
    Bx = EIx.B
    By = EIy.B[:,::-1]
    w_ij_Q1 = np.zeros((len(Bx),len(By)), dtype = complex)
    
    wx = newton_weights(Bx[0], EIx.limits.newtonCotesOrder)
    wy = newton_weights(By[0], EIy.limits.newtonCotesOrder)
    W = np.outer(wx, wy)

    A = A_xy_Q1 * W
    
    for i in range(len(Bx)):
        for j in range(len(By)):
            B_ij_Q1 = np.outer(Bx[i], By[j])
            w_ij_Q1[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q1, A)	
        
            if verbose:
                p.update(n)
                n+=1

    Bx = EIx.B[:,::-1]
    By = EIy.B[:,::-1]
    w_ij_Q2 = np.zeros((len(Bx),len(By)), dtype = complex)
    A = A_xy_Q2 * W
    
    for i in range(len(Bx)):
        for j in range(len(By)):
            B_ij_Q2 = np.outer(Bx[i], By[j])
            w_ij_Q2[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q2, A)
        
            if verbose:
                p.update(n)
                n+=1

    Bx = EIx.B[:,::-1]
    By = EIy.B
    w_ij_Q3 = np.zeros((len(Bx),len(By)), dtype = complex)
    A = A_xy_Q3 * W
    
    for i in range(len(Bx)):
        for j in range(len(By)):
            B_ij_Q3 = np.outer(Bx[i], By[j])
            w_ij_Q3[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q3, A)

            if verbose:
                p.update(n)
                n+=1

    Bx = EIx.B
    By = EIy.B
    w_ij_Q4 = np.zeros((len(Bx),len(By)), dtype = complex)
    A = A_xy_Q4 * W
    
    for i in range(len(Bx)):
        for j in range(len(By)):
            B_ij_Q4 = np.outer(Bx[i], By[j])
            w_ij_Q4[i][j] = dx*dy*np.einsum('ij,ij', B_ij_Q4, A)

            if verbose:
                p.update(n)
                n+=1
                
    if verbose:
        p.finish()
    
    return ROMWeights(w_ij_Q1=w_ij_Q1, w_ij_Q2=w_ij_Q2, w_ij_Q3=w_ij_Q3, w_ij_Q4=w_ij_Q4, EIx=EIx, EIy=EIy)