from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np

class aligo():
	Lambda=1064.0E-9
	etmX_T = 5e-6
	etmX_R = 1.0 - etmX_T
	etmX_r = np.sqrt(etmX_R)
	itmX_T = 0.014
	itmX_R = 1.0 - itmX_T
	itmX_r = np.sqrt(itmX_R)
	etmX_Rc = 2245.0
	itmX_Rc = 1934.0
	LX=3394.5
