import numpy as np

def hermite(n, x): 

	if n == 0: 
        	return 1.0 
	
	elif n == 1:

        	return (2.0 * x) 
	
	elif n == 2:
		
        	return (4.0 * x * x - 2.0) 
	
	elif n == 3:		
            
        	return (8.0 * x * x * x - 12.0 * x) 
	
	
	elif n == 4:
		x_sq = x * x	
        	return (16.0 * x_sq*x_sq - 48.0 * x_sq + 12.0) 
	
	
	elif n == 5:	
		x_sq = x * x 
		x_cb = x * x_sq
        	return (32.0 * x_cb * x_sq - 160.0 * x_cb + 120.0 * x) 
		
	
	elif n == 6:
		x_sq = x * x	
		x_four = x_sq * x_sq
        	return (64.0 * x_four*x_sq - 480.0 * x_four + 720.0 * x_sq - 120.0) 
	
		
	elif n == 7:

		x_sq = x*x
		x_cb = x_sq*x
		x_four = x_cb*x
        	return (128.0 * x_cb*x_four - 1344.0 * x_cb*x_sq + 3360.0 * x_cb - 1680.0 * x) 
	
		
	elif n == 8:
		
		x_sq = x*x
		x_four = x_sq*x_sq
		x_six = x_four*x_sq	
        	return (256.0 * x_four*x_four - 3584.0 * x_six + 13440.0 * x_four - 13440.0 * x_sq + 1680.0) 
	
	elif n == 9:
		x_sq = x*x
		x_cb = x_sq*x
		x_four = x_cb*x	
        	return (512.0 * x_cb*x_cb*x_cb - 9216.0 * x_four*x_cb + 48384.0 * x_cb*x_sq - 80640.0 * x_cb + 30240.0 * x) 
	
	elif n == 10:
	
		x_sq = x*x
		x_cb = x_sq*x
		x_four = x_cb*x
		x_five = x_four * x
	        return (1024.0 * x_five*x_five - 23040.0 * x_four*x_four + 161280.0 * x_cb*x_cb - 403200.0 * x_four + 302400.0 * x_sq - 30240.0) 

	else :

		return (2 * x * hermite(n - 1, x) - 2 * (n - 1) * hermite(n - 2, x))
