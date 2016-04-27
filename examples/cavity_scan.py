# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
"""
Fabry-Perot cavity scan example.

         ]--------(=========================)---->

      laser      ITM       10m cavity      ETM  photodiode

The simulation sets up a parameter list in the form of a Python dictionary,
then populates PyKat with the experimental setup directly.

The cavity is then scanned by tuning the ETM, and the results are plotted.

Note that if you prefer, you can write directly in FINESSE code rather than
using PyKat to build the optical environment - see other examples.

Some terminology:

ITM: initial test mass
ETM: end test mass
HR: highly reflective
AR: anti-reflective

Sean Leavey
s.leavey.1@research.gla.ac.uk

January 2014
"""

import pykat
import pylab as pl
import numpy as np

#######################
# set some parameters #
#######################

parameters = {
	'laser': {
		'power': 30,				# input laser power [W]
		'frequency_offset': 0,
		'phase': 0
	},
	'cavity': {
		'length':	10.8,			# cavity length [m]
		'itm': {				# ITM
			'radius': 0.023225,		# [m]
			'radius_of_curvature': {
				'inner': {		# inner (concave) surface
					'x': -5.7,	# [m]
					'y': -5.7	# [m]
				},
				'outer': {		# outer (convex) surface
					'x': -1.7763,	# [m]
					'y': -1.7763	# [m]
				}
			},
			'thickness': 0.027,		# [m]
			'reflectivity': {		# power reflectivity
				'inner': 0.995,		# inner (concave) surface
				'outer': 0.001		# outer (convex) surface
			},
			'transmission': {		# power transmission
				'inner': 0.005,
				'outer': 0.999
			},
			'tuning_angle': {		# phi
				'inner': 0,
				'outer': 0
			},
			'misalignment': {		# mirror misalignment [rad]
				'inner': {
					'x': 0,
					'y': 0
				},
				'outer': {
					'x': 0,
					'y': 0
				}
			}
		},
		'etm': {				# ETM
			'radius': 0.023225,		# [m]
			'radius_of_curvature': {
				'inner': {		# inner (concave) surface
					'x': 5.7,	# [m]
					'y': 5.7	# [m]
				},
				'outer': {		# outer (convex) surface
					'x': 1.7763,	# [m]
					'y': 1.7763	# [m]
				}
			},
			'thickness': 0.027,		# [m]
			'reflectivity': {		# power reflectivity
				'inner': 0.995,		# inner (concave) surface
				'outer': 0.001		# outer (convex) surface
			},
			'transmission': {		# power transmission
				'inner': 0.005,
				'outer': 0.999
			},
			'tuning_angle': {		# phi
				'inner': 0,
				'outer': 0
			},
			'misalignment': {		# mirror misalignment [rad]
				'inner': {
					'x': 0,
					'y': 0
				},
				'outer': {
					'x': 0,
					'y': 0
				}
			}
		},
	},
	'materials': {
		'bulk': {
			'silica': {
				'refractive_index': 1.45
			}
		}
	}
}

###############################################
# instantiate PyKat object and add components #
###############################################

# instantiate PyKat object
kat = pykat.finesse.kat()

# laser
kat.add(
	pykat.components.laser(
		'laser',					# name
		'n1',						# node
		parameters['laser']['power'],
		parameters['laser']['frequency_offset'],
		parameters['laser']['phase']
	)
)

# add a 1m space between laser and ITM
kat.add(
	pykat.components.space(
		'space1',					# name
		'n1',						# node 1
		'n2',						# node 2
		1						# length [m]
	)
)

##################
# ITM definition #
##################
# This involves three 'components':
#	* a mirror to represent the convex AR surface;
#	* a space representing the thickness of the mirror, with correct refractive index;
#	* a mirror representing the concave HR surface

# AR coating
kat.add(
	pykat.components.mirror(
		'M_ITM_AR',
		'n2',
		'n3',
		parameters['cavity']['itm']['reflectivity']['outer'],
		parameters['cavity']['itm']['transmission']['outer'],
		parameters['cavity']['itm']['tuning_angle']['outer'],
		parameters['cavity']['itm']['radius_of_curvature']['outer']['x'],
		parameters['cavity']['itm']['radius_of_curvature']['outer']['y'],
		parameters['cavity']['itm']['misalignment']['outer']['x'],
		parameters['cavity']['itm']['misalignment']['outer']['y'],
		0,
		parameters['cavity']['itm']['radius'] * 2
	)
)

# bulk mirror material
kat.add(
	pykat.components.space(
		'M_ITM_BULK',
		'n3',
		'n4',
		parameters['cavity']['itm']['thickness'],
		parameters['materials']['bulk']['silica']['refractive_index']
	)
)

# HR coating
kat.add(
	pykat.components.mirror(
		'M_ITM_HR',
		'n4',
		'n5',
		parameters['cavity']['itm']['reflectivity']['inner'],
		parameters['cavity']['itm']['transmission']['inner'],
		parameters['cavity']['itm']['tuning_angle']['inner'],
		parameters['cavity']['itm']['radius_of_curvature']['inner']['x'],
		parameters['cavity']['itm']['radius_of_curvature']['inner']['y'],
		parameters['cavity']['itm']['misalignment']['inner']['x'],
		parameters['cavity']['itm']['misalignment']['inner']['y'],
		0,
		parameters['cavity']['itm']['radius'] * 2
	)
)

##########
# cavity #
##########

kat.add(
	pykat.components.space(
		'space2',
		'n5',
		'n6',
		parameters['cavity']['length']
	)
)

##################
# ETM definition #
##################
# This involves three 'components', just like the ITM definition.

# HR coating
kat.add(
	pykat.components.mirror(
		'M_ETM_HR',
		'n6',
		'n7',
		parameters['cavity']['etm']['reflectivity']['inner'],
		parameters['cavity']['etm']['transmission']['inner'],
		parameters['cavity']['etm']['tuning_angle']['inner'],
		parameters['cavity']['etm']['radius_of_curvature']['inner']['x'],
		parameters['cavity']['etm']['radius_of_curvature']['inner']['y'],
		parameters['cavity']['etm']['misalignment']['inner']['x'],
		parameters['cavity']['etm']['misalignment']['inner']['y'],
		0,
		parameters['cavity']['etm']['radius'] * 2
	)
)

# bulk mirror material
kat.add(
	pykat.components.space(
		'M_ETM_BULK',
		'n7',
		'n8',
		parameters['cavity']['etm']['thickness'],
		parameters['materials']['bulk']['silica']['refractive_index']
	)
)

# AR coating
kat.add(
	pykat.components.mirror(
		'M_ETM_AR',
		'n8',
		'n9',
		parameters['cavity']['etm']['reflectivity']['outer'],
		parameters['cavity']['etm']['transmission']['outer'],
		parameters['cavity']['etm']['tuning_angle']['outer'],
		parameters['cavity']['etm']['radius_of_curvature']['outer']['x'],
		parameters['cavity']['etm']['radius_of_curvature']['outer']['y'],
		parameters['cavity']['etm']['misalignment']['outer']['x'],
		parameters['cavity']['etm']['misalignment']['outer']['y'],
		0,
		parameters['cavity']['etm']['radius'] * 2
	)
)

##############
# photodiode #
##############

# photodiode looking at cavity transmitted light
kat.add(
	pykat.detectors.pd(
		'pd1',
		0,
		'n9'
	)
)

###########################
# Gaussian beam parameter #
###########################

# set q value 1m from ITM, i.e. at the n1 node
# use the utility method for this purpose
kat.space1.n1.q = pykat.optics.gaussian_beams.gauss_param(q = 1.050412 + 24.243836j)
# you can alternatively set w0 and z with gauss_param(w0 = #, z = #)

##############################
# define what we want to see #
##############################

# scan cavity from 0 to 360 degrees
kat.add(pykat.commands.xaxis('lin', [0, 360], kat.M_ETM_HR.phi, 360))

# set maximum TEM mode to model
kat.maxtem = 3

#######################
# run script and plot #
#######################

# run simulation
r = kat.run()

# output the raw FINESSE file that PyKat has generated
scriptList = kat.generateKatScript()
print (''.join(scriptList))

# calculate and print cavity finesse
r1r2 = np.sqrt(parameters['cavity']['itm']['reflectivity']['inner']) * np.sqrt(parameters['cavity']['etm']['reflectivity']['inner'])

finesse = np.pi / (2 * np.arcsin((1 - r1r2) / (2 * np.sqrt(r1r2))))

print ("Cavity finesse: {0:.0f}".format(finesse))

# create plot
pl.plot(r.x, r.y)

# show grid
pl.grid(True)

# set plot limits
pl.xlim((0, 360))

# make plot visible
pl.show()
