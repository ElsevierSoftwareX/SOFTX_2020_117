#import pylab as pl # removed by DDB 1/12/2015
import scipy.optimize as opt
from pykat import finesse
from pykat.detectors import *
from pykat.components import *
from pykat.commands import *
from pykat.structs import *
from numpy import *
# from modematch import modematch
import pykat.optics.ABCD as abcd
import time



def mmf(f, d, q1, q2, c, d2_min):
    '''
    Function used to optimize lenses for modematching
    f - array of slice objects. It gives the allowed focal lengths.
    d - array of slice objects. Gives the allowed range of distances.
    D - total distance between the two beam waists
    q1 - beam parameter at first waist
    q2 - beam parameter at second waist
    '''

    # Optimizes the distances d for each focal lenght in f1 and f2.
    res = opt.brute(mmd, d, args=(f, q1, q2, c, d2_min), \
                    full_output=True, finish=opt.fmin, disp=True)
                    
    # Returns the error of the optimised distances.
    return res[1]


def mmd(d, f, q1, q2, c, d2_min):
    '''
    Function used ot optimize distances between waists and lenses,
    given focal lengths f1 and f2. The parameters are the same as
    defined above for mmf.
    '''
    
    D = d[1]
    d = d[0]
    
    # Checking if the distances are valid 
    if D>d2_min+2*c and d > c and d < D-c-d2_min:
        
        # ABCD-matrices
        M1 = abcd.space(1,d)
        M2 = abcd.lens(f)
        M3 = abcd.space(1,D-d)
        # Total ABCD-matrix
        M = M3*M2*M1

        A = M[0,0]
        B = M[0,1]
        C = M[1,0]
        D = M[1,1]

        # Obtained beam parameter at the second waist position
        q = (A*q1 + B)/(C*q1 + D)
        
        # Returns the difference between the obtained
        # parameter and the target.
        return abs(q-q2)
    
    # If the distances are invalid, return inf
    else:
        return inf


def mmf2(f, d, D, q1, q2, c1, c2, c3):
    '''
    Function used to optimize two lenses in serial for modematching
    f - array of arrays of slice objects giving the allowed focal lengths.
    d - array of slice objects. Gives the allowed range of distances.
    D - total distance between the two beam waists
    q1 - beam parameter at first waist
    q2 - beam parameter at second waist
    '''

    # Focal lengths
    f1 = f[0]
    f2 = f[1]
    # Optimizes the distances d for each focal lenght in f1 and f2.
    res = opt.brute(mmd2, d, args=(f1, f2, D, q1, q2, c1, c2, c3), \
                    full_output=True, finish=opt.fmin)
    # Returns the error of the optimised distances.
    return res[1]


def mmd2(d, f1, f2, L, q1, q2, c1, c2, c3):
    '''
    Function used ot optimize distances between waists and lenses,
    given focal lengths f1 and f2. The parameters are the same as
    defined above for mmf2.
    '''
    
    # Distance waist1 --> lens1
    d1 = d[0]
    # Distance lens1 --> lens2
    d2 = d[1]

    # Checking if the distances are valid 
    if d1 > c1 and d1 < L-c2-c3 and d2 > c2 and \
       d2 < L-c3-d1:

        # ABCD-matrices
        M1 = abcd.space(1,d1)
        M2 = abcd.lens(f1)
        M3 = abcd.space(1,d2)
        M4 = abcd.lens(f2)
        M5 = abcd.space(1,L-d1-d2)
        # Total ABCD-matrix
        M = M5*M4*M3*M2*M1

        A = M[0,0]
        B = M[0,1]
        C = M[1,0]
        D = M[1,1]

        # Obtained beam parameter at the second waist position
        q = (A*q1 + B)/(C*q1 + D)
        
        # Returns the difference between the obtained
        # parameter and the target.
        return abs(q-q2)
    
    # If the distances are invalid, return inf
    else:
        return inf



def moma2(d, f1, f2, D, q1, q2):
    '''
    Unused function that was used to solve a system of two equations
    to modematch. The parameters are the same as defined above for
    mmf and mmd.
    '''
    
    d1 = d[0]
    d2 = d[1]
    
    M1 = abcd.space(1,d1)
    M2 = abcd.lens(f1)
    M3 = abcd.space(1,d2)
    M4 = abcd.lens(f2)
    M5 = abcd.space(1,D-d1-d2)
    
    M = M5*M4*M3*M2*M1
    A = M[0,0]
    B = M[0,1]
    C = M[1,0]
    D = M[1,1]
    
    # q = (A*q1 + B)/(C*q1 + D)

    f = q2 - (A*q1 + B)/(C*q1 + D)
    
    return [f.real, f.imag]



def modematch(q1, q2, f, D, d1, d2_min=.0, c=.01):
    '''
    Simple and not efficient modematching using 1 lens. Quickly translated from the
    two lens case just to test for the filter cavity simulations. Should be fixed
    to a more general, efficient, and useful function at some point.

    q1   d1   f   d2   q2
    |   <-->  |  <-->  |
    | <------ D -----> |

    Input: q1, q2, f, D, d1, d2_min
    q1, q2 - Complex beam parameters to be mode matched.
    f      - Slice-object of available lense. E.g., if the
             available lenses are of focal length 1, 1.5
             and 2 meters, then f = slice(1, 2, 0.5).
    D      - Slice-oject with total distances between q1 and
             q2 that are be searched with brute force. The
             best match of these distances will be optimised
             further. E.g., D = slice(D_min, D_max, D_step).
    d1     - Slice-object with distanes between q1 and the
             lens. E.g., d1 = slice(d1_min, d1_max, d1_step),
             where d1_max < D_max-d2_min-2*c
    d2_min - Minimum distance between the lens and q2
    c      - Minimum distance between the lens and any other
             optic. Here, we are assuming that there are
             optics constraining all the distances, thus c is
             the absolute shortest distance d1 and d2 can take
             even if themselves have 0 as minimum distance.

    Returns: f, d1, d2, D, res
    f      - The focal length of the "best" match found.
    d1     - The d1 of the "best" match found.
    d2     - The d2 of the "best" match found.
    D      - The D of the "best" match found.
    res    - abs(q1-q2) of this match.

    By Daniel Toyra (dtoyra@star.sr.bham.ac.uk)
    '''

    d = [d1, D]
    f = [f]

    # Finding optimal lens
    fsol = opt.brute(mmf, f, args=(d, q1, q2, c, d2_min),
                     full_output=True, finish=None)

    print(fsol)
    f = fsol[0]

    # Finding optimal lens positions given f1, f2
    dsol = opt.brute(mmd, d, args=(f, q1, q2, c, d2_min), \
                     full_output=True, finish=opt.fmin)
    d1 = dsol[0][0]
    D = dsol[0][1]
    d2 = D-d1
    res = dsol[1]

    return f, d1, d2, D, res
    



def modematch2(q1, q2, D, c1, c2, c3):
    '''
    Modematching with 2 lenses. The allowed focal lengths are
    specified inside this function. Later, it would probably
    be better to pass the available lenses as an argument.

    D is total distance between q1 and q2. c1, c2, and c3 are the
    minimum distances between q1 and lens1, lens1 and lens2, and
    lens2 and q2, respectively. c2 also moonlights as the minimum
    distance between any pair of optics. q1 and q2 are the complex
    beam parameters to be matched.

    q1   d1   f1  d2   f2      q2
    |   <-->  |  <-->  |       |
    | <---------- D ---------> |

    By Daniel Toyra (dtoyra@star.sr.bham.ac.uk)
    '''

    print '  Modematching...'
    
    # Controller
    # -----------------------------------------------------
    # Allowed lenses
    f1_range = slice(-1.0, -0.05, 0.05)
    f2_range = slice(0.05, 1.0, 0.05)
    # Number of distance points between each pair of optic
    N = 15
    # -----------------------------------------------------
    
    f = (f1_range, f2_range)
    start1 = c1
    stop1 = D-c2-c3
    step1 = (stop1-start1)/N
    start2 = c2
    stop2 = D-c1-c3
    step2 = (stop2-start2)/N
        
    d1_range = slice(start1, stop1, step1)
    d2_range = slice(start2, stop2, step2)
    d = (d1_range, d2_range)

    # Finding optimal pair of lenses
    fsol = opt.brute(mmf2, f, args=(d, D, q1, q2, c1, c2, c3), \
                     full_output=True, finish=None)
    f1 = fsol[0][0]
    f2 = fsol[0][1]

    # Finding optimal lens positions given f1, f2
    dsol = opt.brute(mmd2, d, args=(f1, f2, D, q1, q2, c1, c2, c3), \
                     full_output=True, finish=opt.fmin)
    d1 = dsol[0][0]
    d2 = dsol[0][1]
    d3 = D-d1-d2
    res = dsol[1]

    # Checking if the solution satisfies the constraints.
    # (it should always be okay, so remove this later if
    # no invalid solutions are found in a while)

   
    if d1 > c1 and d1 < D-c2-c3 and d2 > c2 and \
       d2 < D-d1-c3 and d3 > c3:

        isMM = True

        print '  Match found!'
        print '  res = ' + str(res)
        print '  d1 = ' + str(d1) + ' m'
        print '  d2 = ' + str(d2) + ' m'
        print '  d3 = ' + str(d3) + ' m'
        print '  f1 = ' + str(f1) + ' m'
        print '  f2 = ' + str(f2) + ' m'
        print '  D_tot = ' + str(d1+d2+d3) + ' m'

    else:
        isMM = False
        print '  No match found. Unphysical solution...'
        print d1
        print d2
        print d3
        print res

    
    return f1, f2, d1, d2, d3, res, isMM
    
