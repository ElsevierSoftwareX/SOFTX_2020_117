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
from mpl_toolkits.axes_grid1 import ImageGrid

def main():
	### parameters

	# waist size at laser [m]
	waist = 1e-3

	# power at laser [W]
	power = 1

	# light mode
	mode = 'HG 0 0'

	# grid scale factor
	scale = 2

	# propagation distance [m]
	distance = 1

	### propagation

	# create different grids to demonstrate scaled propagation
	grid1 = oscar.grid(512, 512, 10e-3, 10e-3)
	grid2 = oscar.grid(512, 512, 5e-3, 5e-3)

	# create two identical fields
	field1 = oscar.field(grid1, w=waist, power=power, mode=mode)
	field2 = field1.copy()

	# propagate without scaling
	field1.propagate(distance)

	# propagate with scaling
	field2.scalePropagate(distance, scale, grid2)

	# get magnitudes of the fields
	Z1 = np.abs(field1.amplitude) ** 2
	Z2 = np.abs(field2.amplitude) ** 2

	### plot

	# create figure and image grid for two heatmaps
	fig = pl.figure(figsize=(12, 8))
	imgrid = ImageGrid(fig, 111, nrows_ncols=(1, 2), axes_pad=0.1)

	# figure out global lowest and highest limits of the grids (so we can compare the sizes easily)
	xLim = (min([field1.grid.xaxis.min(), field2.grid.xaxis.min()]), max([field1.grid.xaxis.max(), field2.grid.xaxis.max()]))
	yLim = (min([field1.grid.yaxis.min(), field2.grid.yaxis.min()]), max([field1.grid.yaxis.max(), field2.grid.yaxis.max()]))

	# plot first propagation
	imgrid[0].set_xlim(*xLim)
	imgrid[0].set_ylim(*yLim)
	imgrid[0].set_title('Unscaled Propagation')
	imgrid[0].imshow(Z1, extent=xLim+yLim)

	# plot second propagation
	imgrid[1].set_xlim(*xLim)
	imgrid[1].set_ylim(*yLim)
	imgrid[1].set_title('Scaled Propagation')
	imgrid[1].imshow(Z2, extent=xLim+yLim)

	# show on screen
	pl.show()
if __name__ == '__main__':
	main()

