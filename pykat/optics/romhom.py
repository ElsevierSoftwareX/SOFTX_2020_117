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
from pykat.optics.gaussian_beams import BeamParam
from scipy.linalg import inv
from math import factorial
from pykat.math.hermite import *
from pykat.math import newton_weights
from scipy.integrate import newton_cotes
from multiprocessing import Process, Queue, Array, Value, Event
from pykat.exceptions import BasePyKatException

EmpiricalInterpolant = collections.namedtuple('EmpiricalInterpolant', 'B nodes node_indices limits x worst_error')
ReducedBasis = collections.namedtuple('ReducedBasis', 'RB limits x')
ROMLimits = collections.namedtuple('ROMLimits', 'zmin zmax w0min w0max R mapSamples max_order')

def read_EI(filename, verbose=True):
    import pickle
    
    with open(filename, "rb") as file:
        ei = pickle.load(file)
    
    if verbose:
        print("Map data this ROM was made for in one dimension:")
        print("    Map separation dx = " + str(ei.x[1]-ei.x[0]))
        print("    x range = -{0}m to {0}m".format(max(abs(ei.x))))
        print("    Data points = " + str(ei.x.size * 2))
        print("")
        print("Parameter limits:")
        print("    w0 = {0}m to {1}m".format(ei.limits.w0min, ei.limits.w0max))
        print("    z  = {0}m to {1}m".format(ei.limits.zmin, ei.limits.zmax))
        print("    max order  = {0}".format(ei.limits.max_order))
        print("")
        print("ROM contains {0} basis modes".format(ei.nodes.shape[0]))
        
    return ei
    
                  
class ROMWeights:
    
    def __init__(self, w_ij_Q1, w_ij_Q2, w_ij_Q3, w_ij_Q4, EIx, EIy, nr1=1, nr2=1, direction="reflection_front"):
        self.w_ij_Q1 = w_ij_Q1
        self.w_ij_Q2 = w_ij_Q2
        self.w_ij_Q3 = w_ij_Q3
        self.w_ij_Q4 = w_ij_Q4
        
        self.nr1 = nr1
        self.nr2 = nr2
        
        self.direction = direction 
        
        self.EIx = EIx
        self.EIy = EIy
        
    def writeToFile(self, filename=None, f=None):
        """
        Writes this map's ROM weights to a file
        that can be used with Finesse. The filename
        is appended with '.rom' internally.
        
        Specify either a filename to write the data too, the existing file is overwritten.
        Or provide an open file object to be written too.
        """
        
        if filename is None and f is None:
            raise ValueError("'filename' or open file object 'f' should be specified")
        
        if f is None:
            f = open(filename + ".rom", 'w+')
        
        f.write("direction=%s\n" % self.direction)
        f.write("zmin=%16.16e\n" % self.EIx.limits.zmin)
        f.write("zmax=%16.16e\n" % self.EIx.limits.zmax)
        f.write("w0min=%16.16e\n" % self.EIx.limits.w0min)
        f.write("w0max=%16.16e\n" % self.EIx.limits.w0max)
        f.write("maxorder=%i\n" % self.EIx.limits.max_order)
        f.write("R=%16.16e\n" % self.EIx.limits.R)
        f.write("mapSamples=%i\n" % self.EIx.limits.mapSamples)
        f.write("nr1=%16.16e\n" % self.nr1)
        f.write("nr2=%16.16e\n" % self.nr2)
        
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

        if filename is not None:
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
    if B_matrix is None: return 0
    return np.inner(func[indices].T, B_matrix.T)
   
    
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
    

###################################################################################################
# !!! New ROM code below that doesn't need supercomputer
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


def _write_TS(queue, filename, tssize, Nx, driver):
    from pykat.external.progressbar import ProgressBar
    pb = ProgressBar()
    pb.maxval = tssize
                    
    hfile = h5py.File("%s.h5" % filename, 'a', driver=driver) 
    
    hfile.create_dataset('data', (tssize, Nx), dtype=np.complex128)
    
    i = 0
    
    try:
        while True:
            msg = queue.get()
            
            if msg is None:
                break
            else:
                # Dump each TS into a group
                key = 'TS/%i' % msg[0][4]
                
                hfile[key+"/z"]  = msg[0][0]
                hfile[key+"/w0"] = msg[0][1]
                hfile[key+"/n1"] = msg[0][2]
                hfile[key+"/n2"] = msg[0][3]
                
                hfile["data"][i] = msg[1]
                
                i += 1
                
                if i % 100 == 0:
                    hfile.flush()
                    pb.update(i)
                
    finally:
        hfile.close()
        
def CreateTrainingSetHDF5(filename, maxOrder, z, w0, R, halfMapSamples, NProcesses=1, driver=None):
    
    iq = Queue()
    oq = Queue()
        
    Ns = halfMapSamples
    
    h = R / float(Ns-1) # step size
    
    # Close newton-cotes quadrature goes right upto the boundary
    # unlike previous midpoint rule.
    x = np.linspace(-R, 0, Ns, dtype=np.float64)
    
    w = np.ones(x.shape)
    w[x == 0] = 0.5
    
    nModes = 0
    
    for n in range(0, maxOrder+1):
        for m in range(0, maxOrder+1):
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
    hfile['weights'] = w
    
    hfile.close() # make sure it's closed before
    
    print("Starting processes...")
    # Have a bunch of processes doing the computation and one doing the writing
    iprocesses = [Process(target=_compute_TS, name="irom%i" % i, args=(iq, oq, x, w)) for i in range(NProcesses)]
    oprocess = Process(target=_write_TS, name="orom", args=(oq, filename, tssize, halfMapSamples, driver))
    
    oprocess.start()
    
    try:
        for P in iprocesses:
            P.start()
        
        print("Filling queue...")
        curr = 0
        for n in range(0, maxOrder+1):
            for m in range(0, maxOrder+1):
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
    
    
def _worker_ROM(hdf5Filename, job_queue, result_err, result_idx, event, driver='stdio', chunk=10):
    # h5py drivers: 'core', 'sec2', 'stdio', 'mpio'
    # Need to use something ot her than sec2, the default on OSX,
    # as it doesn't play nice with multiple processes reading files
    
    with h5py.File("%s.h5" % hdf5Filename, driver=driver, mode='r') as file:
        
        while True:
            
            msg = job_queue.get()
            
            if msg is None:
                break
            else:
                TSidx, B, EI_indices = msg
                TSidx = [TSidx[x:x+chunk] for x in range(0, len(TSidx), chunk)]
                
                max_err = 0
                max_idx = -1
                
                for l in TSidx:
                    a = file['data'][l][:]
                    
                    for ll in range(len(a)):
                        res = a[ll] - emp_interp(B, a[ll], EI_indices)
                        _err = np.max( (np.abs(res.real), np.abs(res.imag)) )
                        
                        if _err > max_err or max_idx == -1:
                            max_idx = l[ll]
                            max_err = _err
                    
                result_err.value = max_err
                result_idx.value = max_idx
                
                event.set()

def MakeROMFromHDF5(hdf5Filename, greedyFilename=None, EIFilename=None, tol=1e-10, NProcesses=1, maxRBsize=50, driver="stdio", chunk=10):
    """
    Using a Training Set generated using CreateTrainingSetHDF5 an empirical interpolant is computed.
    
    hdf5Filename = Name of HDF5 file to use
    greedyFilename = Output text file that contains which TS elements were used to make the EI
    EIFilename = Output Pickled file of the EI
    tol = Tolerance for the error on the basis
    NProcesses = Number of processes to use to generate the basis
    maxRBsize = The maximum number of elements in the basis allowed
    driver = The HDF5 driver to use. Use either 'core' or 'stdio'. The former loads the entire TS into memory for each process so can easily overwhelm the computer memory but is faster
    chunk = The number of TS to read from the file at once. Typically 10-100 seem to perform best, faster harddisks can use a larger chunk
    """
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
    EI_indices = []
    RB_matrix = [] 

    V = np.zeros(((maxRBsize), (maxRBsize)), dtype=complex)

    #V[0][0] = RB_matrix[0][EI_indices[0]]
    #invV = inv(V[0:len(EI_indices), 0:len(EI_indices)])
    B = None
    
    RBs = []
    
    if NProcesses > 1:
        queue = Queue()
        locks = [Event() for l in range(NProcesses)]
        
        result_err = [Value('d', np.inf) for l in range(NProcesses)]
        result_idx = [Value('i', -1)     for l in range(NProcesses)]
        
        Names = range(NProcesses)
        procs = [Process(name="process_%i" % l[0], target=_worker_ROM, 
                         args=(hdf5Filename, queue, l[1], l[2], l[3], driver, chunk)) for l in zip(Names, result_err, result_idx, locks)]

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
                       max_order=int(TSdata['maxOrder'].value))
                  
    
    with open(greedyFilename, "w") as f:
        f.write("min w0 = %15.15e\n" % limits.zmin)
        f.write("max w0 = %15.15e\n" % limits.zmax)
        f.write("min z  = %15.15e\n" % limits.w0min)
        f.write("min z  = %15.15e\n" % limits.w0max)
        f.write("R      = %15.15e\n" % limits.R)
        f.write("max order   = %i\n" % limits.max_order)
        f.write("half map samples = %i\n" % limits.mapSamples)
        
        # write initial RB
        _TS = TS[str(next_RB_index)]
        f.write("%15.15e %15.15e %i %i\n" % (_TS["z"].value, _TS["w0"].value, _TS["n1"].value, _TS["n2"].value))
    
        for k in range(1, maxRBsize): 
            _s = time.time()
            
            if NProcesses == 1:
                max_res = []
                TSidx = range(TSsize)
                TSidx = [TSidx[x:x+chunk] for x in range(0, len(TSidx), chunk)]

                for l in TSidx:
                    a = TSdata['data'][l]
                    for ll in range(len(a)):
                        res = a[ll] - emp_interp(B, a[ll], EI_indices)                    
                        max_res.append(np.max( (np.abs(res.real), np.abs(res.imag)) ))

                worst_error = np.max(np.abs(max_res))
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
                print( "Final basis size = %d, Final error = %e, Tolerance=%e" % (k, worst_error, tol) )
                break

            epsilon = TSdata['data'][next_RB_index]
            res = epsilon - emp_interp(B, epsilon, EI_indices)
            
            index_re = np.argmax(abs(res.real))
            index_im = np.argmax(abs(res.imag))
            
            if abs(res.real[index_re]) > abs(res.imag[index_im]):
                index = index_re
            else:
                index = index_im
            
            EI_indices.append(index)
            x_nodes.append(TSdata['x'][index])
            
            print ("worst error = %e at %i on iteration %d" % (worst_error, next_RB_index, k))
            
            RB_matrix.append(res/max(res))
            
            for l in range(len(EI_indices)): # Part of (5) of Algorithm 2: making V_{ij} 
                for m in range(len(EI_indices)): # Part of (5) of Algorithm 2: making V_{ij} 
                    V[m][l] = RB_matrix[l][EI_indices[m]] # Part of (5) of Algorithm 2: making V_{ij}

            invV = inv(V[0:len(EI_indices), 0:len(EI_indices)])
            
            B = B_matrix(invV, np.array(RB_matrix))
            
            _TS = TS[str(next_RB_index)]
            f.write("%15.15e %15.15e %i %i\n" % (_TS["w0"].value, _TS["z"].value, _TS["n1"].value, _TS["n2"].value))
    
            print("Time ", time.time() - _s)
            
    print (time.time() - start, "Seconds")
    
    if NProcesses > 1:
        for P in procs: P.terminate()

    TSdata.close()

    greedyFilenameBase = os.path.splitext(greedyFilename)[0]
    
    print ("Writing to %s" % greedyFilename)
                       
    EI = EmpiricalInterpolant(B=np.matrix(B).real,
                              nodes=np.array(x_nodes).squeeze(),
                              node_indices=np.array(EI_indices).squeeze(),
                              limits=limits,
                              x=x.squeeze(),
                              worst_error=worst_error)
    
    if EIFilename is not None:
        with open("%s.p" % EIFilename, 'wb') as f:
            pickle.dump(EI, f)
        
        print ("Writing to %s.p" % EIFilename)
                              
    return EI
    
    
def makeWeightsNew(smap, EIxFilename, EIyFilename=None, verbose=True, newtonCotesOrderMapWeight=8, direction="reflection_front"):
    
    with open("%s" % EIxFilename, 'rb') as f:
        EIx = pickle.load(f)
        
    if EIyFilename is None:
        EIy = EIx
    else:
        with open("%s" % EIyFilename, 'rb') as f:
            EIy = pickle.load(f)
    
    W_nc = np.outer(newton_weights(smap.x, newtonCotesOrderMapWeight), 
                    newton_weights(smap.y, newtonCotesOrderMapWeight))
                    
    A_xy = smap.z_xy(direction=direction)[::-1, :].T.conj() * W_nc.T
    
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

    wx = np.ones(xm.shape)
    wy = np.ones(xm.shape)
    wx[xm == 0] = 0.5
    wy[ym == 0] = 0.5
    W = np.outer(wx, wy)

    if A_xy_Q1.shape != W.shape or \
            A_xy_Q2.shape != W.shape or \
            A_xy_Q3.shape != W.shape or \
            A_xy_Q4.shape != W.shape:
        raise BasePyKatException("Map data points do not overlap exactly with data points this EI was made for. Consider using the `intepolate=True` option.")

    # make integration weights
    Bx = EIx.B
    By = EIy.B[:,::-1]
    w_ij_Q1 = np.zeros((len(Bx),len(By)), dtype = complex)
    
    A = A_xy_Q1 * W[:,::-1]
    
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
    
    A = A_xy_Q2 * W[::-1,::-1]
    
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

    A = A_xy_Q3 * W[::-1, :]
    
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
    
    return ROMWeights(w_ij_Q1=w_ij_Q1, w_ij_Q2=w_ij_Q2, w_ij_Q3=w_ij_Q3, w_ij_Q4=w_ij_Q4, EIx=EIx, EIy=EIy, direction=direction)
