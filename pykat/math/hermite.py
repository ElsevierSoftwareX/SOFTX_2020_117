from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

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
    elif n == 11:
        return -665280 * x + 2217600 * x**3 - 1774080 * x**5 + 506880 * x**7 - 56320 * x**9 + 2048 * x**11
    elif n == 12:
        return 4096 * x**12-135168 * x**10+1520640 * x**8-7096320 * x**6+13305600 * x**4-7983360 * x**2+665280
    elif n == 13:
        return 8192 * x**13-319488 * x**11+4392960 * x**9-26357760 * x**7+69189120 * x**5-69189120 * x**3+17297280 * x
    elif n == 14:
        return 16384 * x**14-745472 * x**12+12300288 * x**10-92252160 * x**8+322882560 * x**6-484323840 * x**4+242161920 * x**2-17297280
    elif n == 15:
        return 32768 * x**15-1720320 * x**13+33546240 * x**11-307507200 * x**9+1383782400 * x**7-2905943040 * x**5+2421619200 * x**3-518918400 * x
    elif n == 16:
        return 65536 * x**16-3932160 * x**14+89456640 * x**12-984023040 * x**10+5535129600 * x**8-15498362880 * x**6+19372953600 * x**4-8302694400 * x**2+518918400
    elif n == 17:
        return 131072 * x**17-8912896 * x**15+233963520 * x**13-3041525760 * x**11+20910489600 * x**9-75277762560 * x**7+131736084480 * x**5-94097203200 * x**3+17643225600 * x
    elif n == 18:
        return 262144 * x**18-20054016 * x**16+601620480 * x**14-9124577280 * x**12+75277762560 * x**10-338749931520 * x**8+790416506880 * x**6-846874828800 * x**4+317578060800 * x**2-17643225600
    elif n == 19:
        return 524288 * x**19-44826624 * x**17+1524105216 * x**15-26671841280 * x**13+260050452480 * x**11-1430277488640 * x**9+4290832465920 * x**7-6436248698880 * x**5+4022655436800 * x**3-670442572800 * x
    elif n == 20:
        return 1048576 * x**20-99614720 * x**18+3810263040 * x**16-76205260800 * x**14+866834841600 * x**12-5721109954560 * x**10+21454162329600 * x**8-42908324659200 * x**6+40226554368000 * x**4-13408851456000 * x**2+670442572800
    else :
    	return (2 * x * hermite(n - 1, x) - 2 * (n - 1) * hermite(n - 2, x))
