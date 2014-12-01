import numpy as np
from scipy.misc import factorial as fac

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

