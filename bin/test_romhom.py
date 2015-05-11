import pykat
import pylab
import numpy as np
import pickle 
fontsize = 13
labelsize = 12

import matplotlib
matplotlib.rc('xtick', labelsize=labelsize) 
matplotlib.rc('ytick', labelsize=labelsize)
matplotlib.rc('text', usetex=False)

from pykat.optics.maps import surfacemap, mergedmap, aperturemap, read_map
from pykat.optics.romhom import CreateTrainingSetHDF5, MakeROMFromHDF5, makeWeightsNew
from pykat.optics.knm import ROM_HG_knm, square_aperture_HG_knm

from pykat.optics.knm import ROM_HG_knm, newton_weights
from pykat.optics.romhom import u_star_u_mm

from copy import deepcopy

ETM_w0 = (0.0047576005354665598, 0.012037040734172199)
ETM_z = (2110.4274999999998, 2202.60826923077)
ITM_w0 = (0.0047576005354665598, 0.012037040734172199)
ITM_z = (-1884.0725, -1791.89173076923)

# # Make a basis for our problem
# from pykat.optics.romhom import CreateTrainingSetHDF5, MakeROMFromHDF5
#
# Settings for ROM
R = 0.32/2.0
maxOrder = 14
halfMapSamples = 600
tolerance = 1e-13

m = "ETM"

hdf5Filename = m + "_%i_TS" % (maxOrder)
greedyFilename = m + "_%i_greedy" % (maxOrder)
EIFilename = m + "_%i_EI" % (maxOrder)

if m == "ETM":
    z  = np.linspace(ETM_z[0], ETM_z[1], 10)
    w0 = np.linspace(ETM_w0[0], ETM_w0[1], 10)
else:
    z  = np.linspace(ITM_z[0], ITM_z[1], 10)
    w0 = np.linspace(ITM_w0[0], ITM_w0[1], 10)



# CreateTrainingSetHDF5(hdf5Filename, maxOrder, z, w0, R,
#                       halfMapSamples, NProcesses=2)
#
# MakeROMFromHDF5(hdf5Filename, greedyFilename=greedyFilename,
#                 EIFilename=EIFilename, tol=tolerance, NProcesses=4, driver="stdio")
#
# print("DONE!")
# import sys
# sys.exit()








N = 2*halfMapSamples-1
x = np.linspace(-R, R, N)
dx = x[2]-x[1]

mapname = "etm07_s1.map"
a = aperturemap("flat", (N,N), (dx,dx), R)

smap = read_map(mapname)
smap.interpolate(x, x)
smap.data *= 0

mm = mergedmap("flat", (N,N), (halfMapSamples,halfMapSamples), (dx,dx), 1)
#mm.addMap(a)
mm.addMap(smap)

NCs = [0, 1, 5, 10]

w = {}

for NCo in NCs:
    w[NCo] = mm.generateROMWeights("%s.p" % EIFilename, verbose=True, newtonCotesOrder=NCo)
    
qx = pykat.beam_param(w0=w0.mean(), z=z.mean())
qy = qx

mode_i = (0,0)
mode_o = (0,0)

for NCo in NCs:
    fwx = newton_weights(x, NCo)
    
    fw_xy = np.outer(fwx, fwx) * mm.z_xy()
    
    k = square_aperture_HG_knm(mode_i, mode_o, qx, R)
    ROQ = ROM_HG_knm(w[NCo], mode_i, mode_o, qx, qx, qy, qy)
    NC  = dx**2 * np.einsum('ij,ij', fw_xy, np.outer(u_star_u_mm(qx.z, qx.w0, mode_i[0], mode_o[0], x), u_star_u_mm(qy.z, qy.w0, mode_i[1], mode_o[1], x)))
    
    print(NCo, k, ROQ, NC, abs(ROQ-k), abs(NC-k), abs(ROQ-NC))
