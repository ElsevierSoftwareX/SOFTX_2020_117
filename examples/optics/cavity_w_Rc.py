"""
---------------------------------------------------------
Simple example to show the use of the utility functions:
- cavity_w1w2_Rc1Rc2(Lambda, L, w1, w2)
- cavity_info(Lambda, L, Rc1, Rc2)

Andreas Freise 01.12.2014
http://www.gwoptics.org/pykat/
---------------------------------------------------------
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from pykat.optics.fabry_perot import *
import numpy as np


def main():
	# wavelength
	Lambda = 1064.0E-9
	L = 3994.515

	print("Advanced LIGO")
	# attr ITMY Rc -1940.7
	# attr ETMY Rc 2242.4
	w1=0.053
	w2=0.062
	[g1, g2, Rc1, Rc2] = cavity_w1w2_Rc1Rc2(Lambda, L, w1, w2)
	print("w1 = {0}, w2 = {1}".format(w1,w2))
	print("g1 = {0}, Rc1 = {1}".format(g1,Rc1))
	print("g2 = {0}, Rc2 = {1}".format(g2,Rc2))
	print("g1*g2 = {0}".format(g1*g2))

	print("A+")
	# attr ITMX Rc -1852.9
	# attr ETMX Rc 2174.3
	w1=0.08
	w2=0.094
	[g1, g2, Rc1, Rc2] = cavity_w1w2_Rc1Rc2(Lambda, L, w1, w2)
	print("w1 = {0}, w2 = {1}".format(w1,w2))
	print("g1 = {0}, Rc1 = {1}".format(g1,Rc1))
	print("g2 = {0}, Rc2 = {1}".format(g2,Rc2))
	print("g1*g2 = {0}".format(g1*g2))

	# cross check
	[zr,w0,z1,wc1,wc2]=cavity_info(Lambda, L, Rc1, Rc2)
	print("(cross check: w1 = {0}, w2 = {1})".format(wc1,wc2))

		
if __name__ == '__main__':
	main()

