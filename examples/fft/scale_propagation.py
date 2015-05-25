from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pylab as pl

import pykat.oscar as oscar
import numpy as np

def main():
	grid1 = oscar.grid(256, 256, 1e-2, 1e-2, xoffset=0, yoffset=0)
	grid2 = oscar.grid(256, 256, 0.5e-2, 0.5e-2)

	field1 = oscar.field(grid1, w=1e-3, power=1, mode='HG 0 0')
	field2 = oscar.field(grid1, w=1e-3, power=1, mode='HG 0 0')

	field1.propagate(10)
	field2.scalePropagate(10, 1, grid2)

	Z1 = np.abs(field1.amplitude) ** 2
	Z2 = np.abs(field2.amplitude) ** 2

	fig = pl.figure()

	pl.subplot(2, 1, 1)
	ax = pl.gca()

	xLim = (field1.grid.xaxis.min(), field1.grid.xaxis.max())
	yLim = (field1.grid.yaxis.min(), field1.grid.yaxis.max())

	ax.set_xlim(*xLim)
	ax.set_ylim(*yLim)

	im = ax.imshow(Z1, extent=[xLim[0], xLim[1], yLim[0], yLim[1]])
	cb = fig.colorbar(im, ax=ax)

	pl.subplot(2, 1, 2)
	ax = pl.gca()

	xLim = (field2.grid.xaxis.min(), field2.grid.xaxis.max())
	yLim = (field2.grid.yaxis.min(), field2.grid.yaxis.max())

	ax.set_xlim(*xLim)
	ax.set_ylim(*yLim)

	im = ax.imshow(Z2, extent=[xLim[0], xLim[1], yLim[0], yLim[1]])
	cb = fig.colorbar(im, ax=ax)

	pl.show()
if __name__ == '__main__':
	main()

