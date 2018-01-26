import numpy as np

def lens(f):
    return np.matrix([[1,0],[-1.0/f, 1]],dtype=float)

def space(d):
    return np.matrix([[1,d],[0, 1]],dtype=float)

def combine(f1, f2, d=None, verbose = False, q=None): 
    '''
    Computing focal length for the new combined lens.

    The errors are computed assuming that the new lens is placed at the
    location of lens1. The ABCD matrices for propagating the complex beam
    parameter to the location of lens2 are constructed, and the errors are
    computed. 

    Inputs:
    -------
    f1      - Focal length of lens1 [m].
    f2      - Focal length of lens2 [m].
    d       - Distance between the two lenses [m]
    verbose - If true, the result and errors are printed to terminal.
    q       - Complex beam parameter immediately before lens1. Used to compute
              error on the q-value immediately after lens2.

    Returns:
    --------
    f         - Focal length of the compound lens [m].
    rel_errs  - Dictionary containing the relative errors for the elements of the
                ABCD-matrix propgating the beam from f1 to f2, asumming that the
                new thin lens is placed at the location of lens1. If a complex
                beam parameter is given, the dictionary also contains the relative
                errors of propagating this particual beam. 
    '''
    
    if d is None:
        d = 0.0
    # Compound focal length
    f = f1*f2/(f1+f2+d)
    
    # Computing errors
    rel_errs = {}
    if d == 0:
        # No errors if the two thin lenses are located at the same place.
        # Handled like this to avoid divison by zero. 
        rel_errs = np.matrix([[0,0],[0, 0]],dtype=float)
    else:
        # Creates ABCD-matrices for...
        # ...compound lens
        F = lens(f)
        # ...lens1
        F1 = lens(f1)
        # ...lens2
        F2 = lens(f2)
        # ...the space separating the two lenses
        L = space(d)

        # ABCD-matrix for system with 2 lenses
        M = F2*L*F1
        # ABCD-matrix for reduced system with 1 lens located at lens1.
        Mc = L*F
        # Relative errors of ABCD-matrix elements
        rel_errs['ABCD_rdiff'] = np.divide(Mc-M,np.abs(M))
        # Relative errors when propagating the q-value
        if not q is None:
            q2 = (M[0,0]*q + M[0,1]) / (M[1,0]*q + M[1,1])
            q2c = (Mc[0,0]*q + Mc[0,1]) / (Mc[1,0]*q + Mc[1,1])
            rel_errs['q_rdiff'] = np.abs(q2c-q2)/np.abs(q2)
    
    if verbose:
        if not q is None:
            print('q-value relative error: {:.2e}'.format(rel_errs['q_rdiff']) )
        print('ABCD matrix relative errors:')
        print(rel_errs['ABCD_rdiff'])
        
    return f, rel_errs


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('f1', type=float)
    parser.add_argument('f2', type=float)
    parser.add_argument('-d','--d', type = float, required=False)
    parser.add_argument('-q', '--q', type = complex, required=False)
    parser.add_argument('-v', '--verbose', required=False)

    args = vars(parser.parse_args())
    f1 = args['f1']
    f2 = args['f2']
    d = args['d']
    verbose = args['verbose']
    q = args['q']
    if verbose is None:
        verbose = False
        
    f, errs = combine(f1, f2, d=d, verbose = verbose, q=q)
