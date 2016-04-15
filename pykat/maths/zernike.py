import numpy as np
from scipy.misc import factorial as fac
from six.moves import xrange
import math

def zernike_R(m, n, rho):
	
	if ((n-m) % 2):
		return rho*0.0
	
	nnm = (n-m)/2.0
	npm = (n+m)/2.0
	
	R = 0
	
	return sum(rho**(n-2.0*k) * (-1.0)**k * fac(n-k) / ( fac(k) * fac(npm - k) * fac(nnm - k)) for k in xrange(int(nnm) + 1))
	
def zernike(m, n, rho, phi):
	"""
	Computes the zernike polynomial in radial coordinates:
		
		Z_{n}^{m}(rho, phi) = R_{n}^{m}(rho) cos(m * phi) (even)
		                      R_{n}^{m}(rho) sin(m * phi) (odd)
	
	Must satisfy n >= m.
	"""
	
	if (n < 0):
		raise ValueError("n must be larger than 0")
	
	if (abs(m) > n):
		raise ValueError("Must use m <= n")
			
	if m > 0:
		return zernike_R(m, n, rho) * np.cos(m * phi)
	elif m < 0:
		return zernike_R(-m, n, rho) * np.sin(-m * phi)
	else:
		return zernike_R(0, n, rho)


def znm2Rc(A,R):
    '''
    Convertes amplitudes of Zernike polynomials of order n=2 into
    spherical radius of curvature. In case the astigmatic modes
    (m=-1,m=2) are included, the functon returns the maximum and
    minimum curvature.
    
    Inputs: A, R
    A - List of amplitudes of order 2 Zernike polynomials, ordered so that m
        increases with list index. 1 <= len(A) <= 3. [m]
    R - Radius of the mirror surface in the xy-plane. [m]

    Returns: Rc
    Rc - If astigmatic modes are used (len(A) == 2 or 3) a numpy.array of length
         2 containing max and min curvatures is returned. If only the 'defocus'
         mode is used Rc is a float number. [m]

    Based on the Simtools function 'FT_Znm_to_Rc.m' by Charlotte Bond.
    '''
    
    if isinstance(A,list):
        if len(A)==3:
            a20 = A[1]
            a22 = math.sqrt(A[0]**2 + A[2]**2)
            s = np.array([1.0, -1.0])
        elif len(A)==2:
            a20 = 0
            a22 = math.sqrt(A[0]**2 + A[1]**2)
            s = np.array([1.0, -1.0])
        elif len(A)==1:
            a20 = A[0]
            a22 = 0
            s = 0
    elif isinstance(A,float) or isinstance(A,int):
        a20 = A
        a22 = 0
        s = 0
        
    Rc = ((2*a20 + s*a22)**2 + R**2)/(2*(2*a20+s*a22))
    return Rc

def Rc2znm(Rc,R):
    '''
    Converts Radius of curvatue to amplitudes of the second order Zernike
    polynomials.

    Iputs: Rc, R
    Rc    - Radius of curvature. Either a number or a numpy.ndarray with
            minimum and maximum curvature.
    R     - Radius of mirror in xy-plane.

    Returns: A
    A     - Ampltude(s) of the second order Zernike polynomials needed. Is
            a number if Rc is a number, and if R is a numpy.ndarray so is A.
    
    Based on Simtools function 'FT_Rc_to_Znm.m' by Charlotte Bond.
    '''

    # Amplitude in x and y directions
    c = ( Rc - np.sign(Rc)*np.sqrt(Rc**2 - R**2) )/2
    
    if isinstance(Rc, np.ndarray):
        A = np.array([])
        # Adding Z(2,0) amplitude
        A = np.append(A,c.mean())
        # Adding Z(2,2) amplitude
        A = np.append(A,2*(c[0]-A[0]))
    elif isinstance(Rc, float) or isinstance(Rc, int):
        A = c
        
    return A
