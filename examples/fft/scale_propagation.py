from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

"""
Demonstration of oscar.py's scaled propagation method.

Sean Leavey
sean.leavey@ligo.org
May 2015
"""

import pykat.oscar as oscar
import pylab as pl
import numpy as np
#from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib as mpl

def main():
	### parameters

	# power at laser [W]
	power = 1

	# light mode
	mode = 'HG 0 0'

	# waist size of beam at start [m]
	waist = 2e-3

	# unscaled physical grid size [m]
	w0 = 10e-3

	# scaled physical grid size [m]
	w1 = 5e-3

	# propagation distance [m]
	distance = 10

	# grid scale factor
	scale = w1 / w0

	### propagation

	# create different grids to demonstrate scaled propagation
	grid1 = oscar.grid(512, 512, w0, w0)
	grid2 = oscar.grid(512, 512, w1, w1)

	# create input field
	laser = oscar.field(grid1, w=waist, power=power, mode=mode)

	# create three identical fields
	field0 = laser.copy()
	field1 = laser.copy()
	field2 = laser.copy()

	# propagate without scaling
	field1.propagate(distance)

	# propagate with scaling
	field2.scalePropagate(distance, scale, grid2)

	# get magnitudes of the fields
	Z0 = np.abs(field0.amplitude) ** 2
	Z1 = np.abs(field1.amplitude) ** 2
	Z2 = np.abs(field2.amplitude) ** 2

	### plot

	# initial and final physical sizes of grids
	extentInit = [min(field0.grid.xaxis), max(field0.grid.xaxis), min(field0.grid.yaxis), max(field0.grid.yaxis)]
	extentFinal1 = [min(field1.grid.xaxis), max(field1.grid.xaxis), min(field1.grid.yaxis), max(field1.grid.yaxis)]
	extentFinal2 = [min(field2.grid.xaxis), max(field2.grid.xaxis), min(field2.grid.yaxis), max(field2.grid.yaxis)]

	# minimum/maximum values across all signals (for colourmap)
	globalMin = min([Z0.min(), Z1.min(), Z2.min()])
	globalMax = min([Z0.max(), Z1.max(), Z2.max()])

	fig, axes = pl.subplots(nrows=2, ncols=2, figsize=(8, 8))

	# original beams
	axes[0, 0].imshow(Z0, extent=extentInit)
	axes[0, 0].set_title('Original Beam')
	axes[0, 0].set_xlabel('Physical width [m]')
	axes[0, 0].set_ylabel('Physical height [m]')
	axes[0, 1].imshow(Z0, extent=extentInit)
	axes[0, 1].set_title('Original Beam')
	axes[0, 1].set_xlabel('Physical width [m]')
	axes[0, 1].set_ylabel('Physical height [m]')

	# unscaled propagated beam
	axes[1, 0].imshow(Z1, extent=extentFinal1)
	axes[1, 0].set_title('Unscaled Propagated Beam')
	axes[1, 0].set_xlabel('Physical width [m]')
	axes[1, 0].set_ylabel('Physical height [m]')

	# scaled propagated beam
	axes[1, 1].imshow(Z2, extent=extentFinal2)
	axes[1, 1].set_title('Scaled Propagated Beam')
	axes[1, 1].set_xlabel('Physical width [m]')
	axes[1, 1].set_ylabel('Physical height [m]')

	# show on screen
	pl.tight_layout()
	pl.show()
if __name__ == '__main__':
	main()

