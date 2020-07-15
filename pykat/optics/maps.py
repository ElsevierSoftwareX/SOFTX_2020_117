"""
------------------------------------------------------
Utility functions for handling mirror surface
maps. Some functions based on earlier version
in Matlab (http://www.gwoptics.org/simtools/)
Work in progress, currently these functions are
untested!

http://www.gwoptics.org/pykat/
------------------------------------------------------
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from pykat.optics.romhom import makeWeightsNew
from scipy.interpolate import interp2d, interp1d
from scipy.optimize import minimize
from pykat.math.zernike import *        
from pykat.exceptions import BasePyKatException
from copy import deepcopy

import numpy as np
import math
import pickle


class MirrorROQWeights:
    
    def __init__(self, rFront, rBack, tFront, tBack):
        self.rFront = rFront
        self.rBack = rBack
        self.tFront = tFront
        self.tBack = tBack
    
    def writeToFile(self, romfilename):
        with open(romfilename + ".rom", "w+") as f:
            if self.rFront is not None: self.rFront.writeToFile(f=f)
            if self.rBack  is not None: self.rBack.writeToFile(f=f)
            if self.tFront is not None: self.tFront.writeToFile(f=f)
            if self.tBack  is not None: self.tBack.writeToFile(f=f)
            
class surfacemap(object):
    
    def __init__(self, name, maptype, size=None, center=None, step_size=(1,1), scaling=1.0e-9, data=None,
                 notNan=None, zOffset=None, xyOffset=(.0,.0)):
        '''
        size, center, step_size, xyOffset are all tuples of the form (x, y),
        i.e., (col, row).
        '''
        assert(len(step_size) == 2)
        assert(step_size[0] > 0 and step_size[1] > 0)
        assert(len(xyOffset) == 2)
        assert(size is not None and len(size) == 2)
        assert(center is not None and len(center) == 2)
        
        self.name = name
        self.type = maptype
        
        if data is None:
            if size is None:
                raise BasePyKatException("One of the parameters data or size needs to be specified to create map")
            elif isinstance(size, tuple) or isinstance(size, list):
                if len(size) == 2:
                    self.data = np.zeros(size[::-1])
                else:
                    raise BasePyKatException("Parameter size must have length 2")
            elif isinstance(size, int):
                self.data = np.zeros((size, size))
        else:
            self.data = data

        self.notNan = notNan

        # Currently "beam center", i.e., mirror_center + xyOffset. 
        self.step_size = step_size
        self.center = center
        self.scaling = scaling
        self._RcRemoved = None
        
        # Offset of fitted sphere. Proably unnecessary to have here.
        self.zOffset = zOffset
        self.__interp = None
        self._zernikeRemoved = {}
        self._betaRemoved = None
        self._xyOffset = xyOffset

        self._rom_weights = None
        
    def write_map(self, filename, binary=False):
        with open(filename,'w') as mapfile:
            
            mapfile.write("% Surface map\n")
            mapfile.write("% Name: {0}\n".format(self.name))
            mapfile.write("% Type: {0}\n".format(self.type))
            mapfile.write("% Size: {0} {1}\n".format(self.data.shape[0], self.data.shape[1]))
            mapfile.write("% Optical center (x,y): {0} {1}\n".format(self.center[0], self.center[1]))
            mapfile.write("% Step size (x,y): {0} {1}\n".format(self.step_size[0], self.step_size[1]))
            mapfile.write("% Scaling: {0}\n".format(float(self.scaling)))
            mapfile.write("\n\n")

            for i in range(0, self.data.shape[0]):
                for j in range(0, self.data.shape[1]):
                    mapfile.write("%.15g " % self.data[i,j])
                mapfile.write("\n")

    @property
    def xyOffset(self):
        '''
        The offset from the mirror center measured in meters.
        '''
        return self._xyOffset

    @xyOffset.setter
    def xyOffset(self,offset):
        '''
        Give 'offset' is in meters.
        '''
        
        # These two lines aren't compatible with recenter(). Changed.
        # x0new = self.center[0] + (-self.xyOffset[0] + offset[0])/self.step_size[0]
        # y0new = self.center[1] + (-self.xyOffset[1] + offset[1])/self.step_size[1]
        
        x0new = self.center[0] + offset[0]/self.step_size[0]
        y0new = self.center[1] + offset[1]/self.step_size[1]
        
        # Checking so new suggested "center" is within the the mirror surface
        if (x0new>0 and x0new<self.size[0]-1 and y0new>0 and y0new<self.size[1]-1 and
            self.notNan[int(round(y0new)), int(round(x0new))]):
            
            self.center = (x0new, y0new)
            self._xyOffset = (offset[0], offset[1])
        
        else:
            raise BasePyKatException( ('Error in xyOffset: ({:.2e}, {:.2e}) m --> pos = ({:.2f}, {:.2f}) '+
                  'is outside mirror surface.').format(offset[0],offset[1],x0new,y0new))
                    
    @property
    def betaRemoved(self):
        return self._betaRemoved
    
    @betaRemoved.setter
    def betaRemoved(self, v):
        self._betaRemoved = v
        
    @property
    def zernikeRemoved(self):
        return self._zernikeRemoved
    
    @zernikeRemoved.setter
    def zernikeRemoved(self, v):
        '''
        v = tuple(m,n,amplitude)
        '''
        self._zernikeRemoved["%i%i" % (v[0],v[1])] = v

    @property
    def data(self):
        return self.__data
    
    @data.setter
    def data(self, value):
        self.__data = value
        self.__interp = None

    @property
    def notNan(self):
        return self.__notNan

    @notNan.setter
    def notNan(self,value):
        if value is None:
            if not hasattr(self,'notNan'):
                self.__notNan = np.ones(self.size[::-1], dtype=bool)
        else:
            self.__notNan = value
            
    @property
    def center(self):
        return self.__center
    
    @center.setter
    def center(self, value):
        if value is None:
            if not hasattr(self,'center'):
                self.__center = self.recenter()
        elif isinstance(value, tuple):
            if len(value) == 2:
                self.__center = value
        elif isinstance(value, list):
            if len(value) == 2:
                self.__center = (value[0],value[1])
        elif isinstance(value, np.ndarray):
            if len(value) == 2:
                self.__center = (value[0],value[1])
        elif isinstance(value, float) or isinstance(value, int):
            self.__center = (value, value)
        else:
            raise BasePyKatException("Invalid format of center, array of length 2 wanted.")
            
        self.__interp = None
    
    @property
    def step_size(self):
        return self.__step_size
    
    @step_size.setter
    def step_size(self, value):
        if isinstance(value, tuple):
            self.__step_size = value
        elif isinstance(value, list):
            self.__step_size = (value[0],value[1])
        elif isinstance(value, float) or isinstance(value, int):
            self.__step_size = (value,value)
        self.__interp = None

    @property
    def scaling(self):
        return self.__scaling
    
    @scaling.setter
    def scaling(self, value):
        self.__scaling = value
        self.__interp = None

    # CHANGED! I swapped: self.data.shape[0]  <----> self.data.shape[1]. Because
    # the rows/columns in the data matrix are supposed to be y/x-values(?).
    # This may mess things up in other methods though... I fixed the plot-method.
    # /DT
    @property
    def x(self):
        return self.step_size[0] * (np.array(range(0, self.data.shape[1])) - self.center[0])
    
    @property
    def y(self):
        return self.step_size[1] * (np.array(range(0, self.data.shape[0]))- self.center[1])
        
    # CHANGED! Since everything else (step_size, center) are given as (x,y), and not
    # as (row, column), I changed this to the same format. /DT
    @property
    def size(self):
        return self.data.shape[::-1]
        
    @property
    def offset(self):
        return np.array(self.step_size)*(np.array(self.center) - (np.array(self.size)-1)/2.0)
    
    @property
    def ROMWeights(self):
        return self._rom_weights
    
    def z_xy(self, x=None, y=None, wavelength=1064e-9, direction="reflection_front", nr1=1.0, nr2=1.0):
        """
        For this given map the field perturbation is computed. This data
        is used in computing the coupling coefficient. It returns a grid
        of complex values representing the change in amplitude or phase
        of the field.
        
            x, y      : Points to interpolate at, 'None' for no interpolation.
            
            wavelength: Wavelength of light in vacuum [m]
            
            direction : Sets which distortion to return, as beams travelling
                        in different directions will see different distortions.
                        Options are:
                                "reflection_front"
                                "transmission_front" (front to back)
                                "transmission_back" (back to front)
                                "reflection_back"
                                
            nr1       : refractive index on front side
            
            nr2       : refractive index on back side
            
        """
        
        assert(nr1 >= 1)
        assert(nr2 >= 1)
        
        if x is None and y is None:
            data = self.scaling * self.data
        else:
            if self.__interp is None:
                self.__interp = interp2d(self.x, self.y, self.data * self.scaling)
                
            data = self.__interp(x, y)
        
        if direction == "reflection_front" or direction == "reflection_back":
            if "phase" in self.type:
                k = math.pi * 2 / wavelength
                
                if direction == "reflection_front":
                    return np.exp(-2j * nr1 * k * data)
                else:
                    return np.exp(2j * nr2 * k * data[:,::-1])
                
            elif "absorption" in self.type:
                if direction == "reflection_front":
                    return np.sqrt(1.0 - data)
                else:
                    return np.sqrt(1.0 - data[:, ::-1])
            elif "surface_motion" in self.type:
                if direction == "reflection_front":
                    return data
                else:
                    return data[:, ::-1]
            else:
                raise BasePyKatException("Map type needs handling")
                
        elif direction == "transmission_front" or direction == "transmission_back":
            if "phase" in self.type:
                k = math.pi * 2 / wavelength
                
                if direction == "transmission_front":
                    return np.exp((nr1-nr2) * k * data)
                else:
                    return np.exp((nr2-nr1) * k * data[:, ::-1])
                
            elif "absorption" in self.type:
                if direction == "transmission_front":
                    return np.sqrt(1.0 - data)
                else:
                    return np.sqrt(1.0 - data[:, ::-1])
                    
            elif "surface_motion" in self.type:
                return np.ones(data.shape)
            else:
                raise BasePyKatException("Map type needs handling")
                
        else:
            raise ValueError("Direction not valid")
        

    
    def generateROMWeights(self, EIxFilename, EIyFilename=None, nr1=1.0, nr2=1.0, verbose=False, interpolate=False, newtonCotesOrder=8):
        
        if interpolate == True:
            # Use EI nodes to interpolate if we
            with open(EIxFilename, 'rb') as f:
                EIx = pickle.load(f)

            if EIyFilename is None:
                EIy = EIx
            else:
                with open(EIyFilename, 'rb') as f:
                    EIy = pickle.load(f)

            x = EIx.x
            x.sort()
            nx = np.unique(np.hstack((x, -x[::-1])))
        
            y = EIy.x
            y.sort()
            ny = np.unique(np.hstack((y, -y[::-1])))
            
            self.interpolate(nx, ny)
        
        w_refl_front, w_refl_back, w_tran_front, w_tran_back = (None, None, None, None)
        
        if "reflection" in self.type or "both" in self.type:
            w_refl_front = makeWeightsNew(self, EIxFilename, EIyFilename,
                                      verbose=verbose, newtonCotesOrderMapWeight=newtonCotesOrder,
                                      direction="reflection_front")
            
            w_refl_front.nr1 = nr1
            w_refl_front.nr2 = nr2
            
            w_refl_back = makeWeightsNew(self, EIxFilename, EIyFilename,
                                      verbose=verbose, newtonCotesOrderMapWeight=newtonCotesOrder,
                                      direction="reflection_back")
            
            w_refl_back.nr1 = nr1
            w_refl_back.nr2 = nr2

        if "transmission" in self.type or "both" in self.type:                                      
            w_tran_front = makeWeightsNew(self, EIxFilename, EIyFilename,
                                      verbose=verbose, newtonCotesOrderMapWeight=newtonCotesOrder,
                                      direction="transmission_front")

            w_refl_front.nr1 = nr1
            w_refl_front.nr2 = nr2
                                            
            w_tran_back  = makeWeightsNew(self, EIxFilename, EIyFilename,
                                      verbose=verbose, newtonCotesOrderMapWeight=newtonCotesOrder,
                                      direction="transmission_back")
            
            w_refl_back.nr1 = nr1
            w_refl_back.nr2 = nr2
            
        self._rom_weights = MirrorROQWeights(w_refl_front, w_refl_back, w_tran_front, w_tran_back)
        
        return self._rom_weights
            
    def interpolate(self, nx, ny, **kwargs):
        """
        Interpolates the map for some new x and y values.
        
        Uses scipy.interpolate.interp2d and any keywords arguments are
        passed on to it, thus settings like interpolation type and
        fill values can be set.
        
        The range of nx and ny must contain the value zero so that the
        center point of the map can be set.
        """

        D = interp2d(self.x, self.y, self.data, **kwargs)
        
        data = D(nx-self.offset[0], ny-self.offset[1])
        
        Dx = interp1d(nx, np.arange(0,len(nx)))
        Dy = interp1d(ny, np.arange(0,len(ny)))
        
        self.center = (Dx(0), Dy(0))
        self.step_size = (nx[1]-nx[0], ny[1]-ny[0])
        self.data = data

    # xlim and ylim given in centimeters
    def plot(self, show=True, clabel=None, xlim=None, ylim=None, isBlock=False):

        import matplotlib
        import matplotlib.pyplot as plt
        
        if xlim is not None:
            # Sorts out the x-values within xlim
            _x = np.logical_and(self.x<=max(xlim)/100.0, self.x>=min(xlim)/100.0)
            xmin = np.min(np.where(_x == True))
            xmax = np.max(np.where(_x == True))
        else:
            # Uses the whole available x-range
            xmin = 0
            xmax = len(self.x)-1
            xlim = [self.x.min()*100, self.x.max()*100]
    
        if ylim is not None:
            # Sorts out the y-values within ylim
            _y = np.logical_and(self.y<=max(ylim)/100.0, self.y>=min(ylim)/100.0)
            ymin = np.min(np.where(_y == True))
            ymax = np.max(np.where(_y == True))
        else:
            # Uses the whole available y-range
            ymin = 0
            ymax = len(self.y)-1
            ylim = [self.y.min()*100, self.y.max()*100]
            
        # CHANGED! y corresponds to rows and x to columns, so this should be
        # correct now. Switch the commented/uncommented things to revert. /DT
        # ------------------------------------------------------
        # min and max of z-values
        zmin = self.data[ymin:ymax,xmin:xmax].min()
        zmax = self.data[ymin:ymax,xmin:xmax].max()
        # zmin = self.data[xmin:xmax,ymin:ymax].min()
        # zmax = self.data[xmin:xmax,ymin:ymax].max()
        # ------------------------------------------------------
        
        # 100 factor for scaling to cm
        xRange = 100*self.x
        yRange = 100*self.y
        # xRange, yRange = np.meshgrid(xRange,yRange)
        
        fig = plt.figure()

        # Important to remember here is that xRange corresponds to column and
        # yRange to row indices of the matrix self.data.
        pcm = plt.pcolormesh(xRange, yRange, self.data)
        pcm.set_rasterized(True)
        
        plt.xlabel('x [cm]')
        plt.ylabel('y [cm]')

        if xlim is not None: plt.xlim(xlim)
        if ylim is not None: plt.ylim(ylim)

        plt.title('Surface map {0}, type {1}'.format(self.name, self.type))
        
        cbar = plt.colorbar()
        cbar.set_clim(zmin, zmax)
        # cbar.set_clim(-1.86, 1.04)
        
        if clabel is not None:
            cbar.set_label(clabel)
    
        if show:
            plt.show(block=isBlock)
        
        return fig

    def rescaleData(self,scaling):
        '''
        Rescaling data. Assumes that the current map is scaled correctly, i.e.,
        surfacemap.scaling*surfacemap.data would give the data in meters.

        Input: scaling
        scaling  - float number that representing the new scaling. E.g.,
                   scaling = 1.0e-9 means that surfacemap.data is converted to
                   being given in nano meters.
        '''
        if isinstance(scaling, float):
            self.data[self.notNan] = self.data[self.notNan]*self.scaling/scaling
            self.scaling = scaling
        else:
            print('Error: Scaling factor must be of type float, {:s} found'.
                  format(type(scaling)))
            

    def rms(self, w=None):
        '''
        Returns the root mean square value of the mirror map, expressed in the unit
        given by self.scaling. If w [m] is specified, only the area inside the radius
        is considered. The area is r<w is centered around the self.center, i.e., the
        beam center (not mirror center).
        '''
        if w is None:
            return math.sqrt((self.data[self.notNan]**2).sum()/self.notNan.sum())
        else:
            R = self.find_radius(unit='meters')
            if w>=R:
                return math.sqrt((self.data[self.notNan]**2).sum()/self.notNan.sum())
            else:
                rho = self.createPolarGrid()[0]
                inside = np.zeros(self.data.shape,dtype=bool)
                tmp = rho<w
                inside[tmp] = self.notNan[tmp]
                return math.sqrt((self.data[inside]**2).sum()/inside.sum())
            
    def avg(self, w=None):
        '''
        Returns the average height of the mirror map, expressed in the unit
        given by self.scaling. If w [m] is specified, only the area inside the radius
        is considered. The area is r<w is centered around the self.center, i.e., the
        beam center (not mirror center).
        '''
        if w is None:
            tot = self.data[self.notNan].sum()
            sgn = np.sign(tot)
            return sgn*math.sqrt(sgn*tot)/self.notNan.sum()
        else:
            R = self.find_radius(unit='meters')
            if w>=R:
                tot = self.data[self.notNan].sum()
                sgn = np.sign(tot)
                return sgn*math.sqrt(sgn*tot)/self.notNan.sum()
            else:
                rho = self.createPolarGrid()[0]
                inside = np.zeros(self.data.shape,dtype=bool)
                tmp = rho<w
                inside[tmp] = self.notNan[tmp]
                tot = self.data[inside].sum()
                sgn = np.sign(tot)
                return sgn*math.sqrt(sgn*tot)/inside.sum()
            
            
    def find_radius(self, method = 'max', unit='points'):
        '''
        Estimates the radius of the mirror in the xy-plane. If xyOffset is different from
        (0,0), then the radius will be the minimum distance from the beam center to the
        nearest mirror edge.
        
        method   - 'max' gives maximal distance from beam centre to the edge.
                 - 'xy' returns the distance from centre to edge along x- and y-directions.
                   Doesn't make sense if the beam center is offset from the mirror center.
                 - 'area' calculates the area and uses this to estimate the mean radius.
                 - 'min' gives minimal distance from beam center to the edge
        unit     - 'points' gives radius in data points.
                   'meters' gives radius in meters.

        Based on Simtools function 'FT_find_map_radius.m' by Charlotte Bond.
        '''

        # Row and column indices of non-NaN
        rIndx, cIndx = self.notNan.nonzero()
        if unit == 'meters' or unit == 'metres' or unit=='Meters' or unit=='Metres':
            # Offsets so that (0,0) is in the center.
            x = self.step_size[0]*(cIndx - self.center[0])
            y = self.step_size[1]*(rIndx - self.center[1])
        else:
            # Offsets so that (0,0) is in the center.
            x = cIndx - self.center[0]
            y = rIndx - self.center[1]
            
        # Maximum distance from center to the edge of the mirror.
        if method=='max' or method=='Max' or method=='MAX':
            r = math.sqrt((x**2 + y**2).max())
        # x and y radii
        elif method=='xy' or method=='XY' or method == 'yx' or method == 'YX':
            r = []
            r.append(abs(x).max()+0.5)
            r.append(abs(y).max()+0.5)
        # Mean radius by dividing the area by pi. Should change this in case
        # x and y step sizes are different.
        elif method=='area' or method=='Area' or method=='AREA':
            if unit == 'meters' or unit == 'metres' or unit=='Meters' or unit=='Metres':
                r = self.step_size[0]*math.sqrt(len(cIndx)/math.pi)
            else:
                r = math.sqrt(len(cIndx)/math.pi)
        elif method=='min':
            # Arrays of row and column indices where data is NaN.
            r,c = np.where(self.notNan == False)
            # Sets the map center to index [0,0]
            x = c-self.center[0]
            y = r-self.center[1]
            if unit == 'meters' or unit == 'metres' or unit=='Meters' or unit=='Metres':
                x = x*self.step_size[0]
                y = y*self.step_size[1]
            # Minimum radius squared
            r = math.sqrt((x**2+y**2).min())
        return r


    def zernikeConvol(self,n,m=None):
        '''
        Performs a convolution with the map surface and Zernike polynomials.

        Input: n, m
        n  -  If m (see below) is set to None, n is the maximum order of Zernike
              polynomials that will be convolved with the map. If m is an integer
              or a list, only the exact order n will be convolved.
        m  -  = None. All Zernike polynomials up to and equal to order n will
              be convolved.
              = value [int]. The Zernike polynomial characterised by (n,m) will
              be convolved witht the map.
              = list of integers. The Zernike polynomials of order n with m-values
              in m will be convolved with the map.
             
        Returns: A
        A  -  If m is None:
                List with lists of amplitude coefficients. E.g. A[n] is a list of
                amplitudes for the n:th order with m ranging from -n to n, i.e.,
                A[n][0] is the amplitude of the Zernike polynomial Z(m,n) = Z(-n,n).
              If m is list or int:
                Returns list with with Zernike amplitudes ordered in the same way as
                in the input list m. If m == integer, there is only one value in the
                list A.
        
        Based on the simtool function 'FT_zernike_map_convolution.m' by Charlotte
        Bond.
        '''
        
        mVals = []
        if m is None:
            nVals = range(0,n+1)
            for k in nVals:
                mVals.append(range(-k,k+1,2))
        elif isinstance(m, list):
            nVals = range(n,n+1)
            mVals.append(m)
        elif isinstance(m,int):
            nVals = range(n,n+1)
            mVals.append(range(m,m+1))
        else:
            raise Exception("Unhandled `m` argument: %s" % m)
            
        # Amplitudes
        A = []
        # Radius
        R = self.find_radius(unit='meters')
        
        # Grid for for Zernike-polynomials (same size as map).
        rho, phi = self.createPolarGrid()
        rho = rho/R
        # Loops through all n-values up to n_max
        for k in range(len(nVals)):
            n = nVals[k]
            # Loops through all possible m-values for each n.
            tmp = []
            for m in mVals[k]:
                # (n,m)-Zernike polynomial for this grid.
                Z = zernike(m, n, rho, phi)
                # Z = znm(n, m, rho, phi)
                # Convolution (unnecessary conjugate? map.data is real numbers.)
                c = (Z[self.notNan]*np.conjugate(self.data[self.notNan])).sum()
                cNorm = (Z[self.notNan]*np.conjugate(Z[self.notNan])).sum()
                c = c/cNorm
                tmp.append(c)
            A.append(tmp)
        if len(A) == 1:
            A = A[0]
            if len(A)==1:
                A = A[0]
                
        return A

    def rmZernikeCurvs(self, zModes='all',interpol=False):
        '''
        Removes curvature from mirror map by fitting Zernike polynomials to the 
        mirror surface. Also stores the estimated radius of curvature
        and the Zernike polynomials and amplitudes.
        
        zModes   - 'defocus' uses only Zernike polynomial (n,m) = (2,0), which has a
                   parabolic shape.
                   'astigmatism' uses Zernike polynomials (n,m) = (2,-2) and (n,m) = (2,2).
                   'all' uses both defocus and the astigmatic modes.
        interpol - Boolean where 'true' means that the surface map is interpolated to get
                   more data points.

        Returns [Rc,A]
        Rc       - Estimated radius of curvature removed [m].
        A        - Amplitudes of Zernike polynomials removed [surfacemap.scaling].

        Based on the Simtools functions 'FT_remove_zernike_curvatures_from_map.m',
        'FT_zernike_map_convolution.m', etc. by Charlotte Bond.
        '''
        
        tmp = deepcopy(self)
        ny,nx = self.data.shape
        # Interpolating for more precise convolution, in case size of map is small. Not sure
        # this is necessary. It isn't used in the latest Simtools tutorials.
        if interpol and (nx<1500 or ny<1500):
            # Number of extra steps inserted between two adjacent data points.
            N = math.ceil(1500.0/min(nx,ny))
            # New arrays of x-values and y-values
            x = np.arange(tmp.x[0],tmp.x[-1]+tmp.step_size[0]/(N+1),tmp.step_size[0]/N)
            y = np.arange(tmp.y[0],tmp.y[-1]+tmp.step_size[1]/(N+1),tmp.step_size[1]/N)

            '''
            # --------------------------------------------------------------
            # Interpolating data
            tmp.interpolate(x,y)
            tmp.plot()
            # Interpolating for notNan.
            g = interp2d(self.x, self.y, self.notNan, kind='linear')
            tmp.notNan = np.round(g(x,y)) # round() or floor() here?! Can't decide...
            # Converting binary to boolean
            tmp.notNan = tmp.notNan==1
            # --------------------------------------------------------------
            '''
            # Later when compatible with interpolation function, switch code
            # below for code above (but need to change 
            # --------------------------------------------------------------
            # Interpolation functions, for the data (f) and for notNan (g).
            f = interp2d(tmp.x, tmp.y, tmp.data, kind='linear')
            g = interp2d(tmp.x, tmp.y, tmp.notNan, kind='linear')
            # Interpolating
            tmp.data = f(x,y)

            tmp.notNan = np.round(g(x,y)) # round() or floor() here?! Can't decide...
            # Converting binary to boolean
            tmp.notNan = tmp.notNan==1

            # Setting new step size
            tmp.step_size = (x[1]-x[0],y[1]-y[0])
            # Setting new center
            fx = interp1d(x, np.arange(0,len(x)))
            fy = interp1d(y, np.arange(0,len(y)))
            tmp.center = (fx(0.0), fy(0.0))
            # --------------------------------------------------------------

        # Convolution between the surface map and Zernike polynomials
        # to get amplitudes of the polynomials. Only interested in n=2.
        A = tmp.zernikeConvol(2)[2]
        
        # Radius of mirror.
        R = self.find_radius(unit='meters')
        # Grid in polar coordinates
        rho,phi = self.createPolarGrid()
        rho = rho/R
        
        # Creates the choosen Zernike polynomials and removes them from the
        # mirror map.
        if zModes=='all' or zModes=='All':
            for k in range(0,3):
                m = (k-1)*2
                Z = A[k]*zernike(m, 2, rho, phi)
                self.data[self.notNan] = self.data[self.notNan]-Z[self.notNan]
                self.zernikeRemoved = (m, 2, A[k])
            # Estimating radius of curvature
            Rc = znm2Rc([a*self.scaling for a in A], R)
        elif zModes=='astigmatism' or zModes=='Astigmatism':
            for k in range(0,3,2):
                m = (k-1)*2
                Z = A[k]*zernike(m, 2, rho, phi)
                self.data[self.notNan] = self.data[self.notNan]-Z[self.notNan]
                self.zernikeRemoved = (m, 2, A[k])
            Rc = znm2Rc([a*self.scaling for a in A[::2]], R)
        elif zModes=='defocus' or zModes=='Defocus':
            Z = A[1]*zernike(0, 2, rho, phi)
            self.data[self.notNan] = self.data[self.notNan]-Z[self.notNan]
            self.zernikeRemoved = (0, 2, A[1])
            Rc = znm2Rc(A[1]*self.scaling, R)
        
        self.RcRemoved = Rc
        
        return self.RcRemoved, self.zernikeRemoved


    def rmTilt(self, method='fitSurf', w=None, xbeta=None, ybeta=None, zOff=None):
        
        R = self.find_radius(unit='meters')
        
        if method=='zernike' or method=='Zernike':
            A1 = self.rmZernike(1, m = 'all')
            # Calculating equivalent tilts in radians
            xbeta = np.arctan(A1[1]*self.scaling/R)
            ybeta = np.arctan(A1[0]*self.scaling/R)
            self.betaRemoved = (xbeta,ybeta)
            return A1, xbeta, ybeta

        else:
            X,Y = np.meshgrid(self.x,self.y)
            r2 = X**2 + Y**2

            if w is not None:
                weight = 2/(math.pi*w**2) * np.exp(-2*r2[self.notNan]/w**2)

            def f(p):
                # This is used in simtools, but gives less good results from what I have seen, also in simtools.
                # I would guess the problem is that the initial guess is 0, and the minimising algorithm would
                # need to go quite far away from zero to find a minumum even slightly away way from zero.
                # I think the idea is to make the algorithm to search for small angles. Of course, if used, this
                # must be converted back after minimisation.
                # p[0] = p[0]*1.0e-9
                # p[1] = p[1]*1.0e-9
                Z = self.createSurface(0,X,Y,p[2],0,0,p[0],p[1])
                if w is None:
                    res = math.sqrt(((self.data[self.notNan]-Z[self.notNan])**2).sum())/self.notNan.sum()
                else:
                    # weight = 2/(math.pi*w**2) * np.exp(-2*r2[self.notNan]/w**2)
                    res = math.sqrt((weight*(self.data[self.notNan]-Z[self.notNan])**2).sum())/weight.sum()

                return res

            if xbeta is None:
                xbeta = 0
            if ybeta is None:
                ybeta = 0
            if zOff is None:
                zOff = 0

            params = [xbeta,ybeta,zOff]

            opts = {'xtol': 1.0e-9, 'ftol': 1.0e-9, 'maxiter': 1000, 'maxfev': 1000, 'disp': False}
            out = minimize(f, params, method='Nelder-Mead', options=opts)
            if not out['success']:
                msg = '  Warning: ' + out['message'].split('.')[0] + ' (nfev={:d}).'.format(out['nfev'])
                print(msg)

            xbeta = out['x'][0]
            ybeta = out['x'][1]
            zOff = out['x'][2]

            Z = self.createSurface(0,X,Y,zOff,0,0,xbeta,ybeta)
            self.data[self.notNan] = self.data[self.notNan] - Z[self.notNan]
            self.betaRemoved = (xbeta, ybeta)

            # Equivalent Zernike amplitude
            A1 = R*np.tan(np.array([ybeta,xbeta]))/self.scaling
            self.zernikeRemoved = (-1,1,A1[0])
            self.zernikeRemoved = (1,1,A1[1])
            
            return A1,xbeta,ybeta,zOff
        

    def rmSphericalSurf(self, Rc0, w=None, zOff=None, isCenter=[False,False], maxfev=2000):
        '''
        Fits spherical surface to the mirror map and removes it.

        Inputs: Rc0, w, zOff, isCenter
        Rc       - Initial guess of the Radius of curvature. [m]
        w        - Gaussian weighting parameter. The distance from the center where the
                   weigths have decreased by a factor of exp(-2) compared to the center.
                   Should preferrably be equal to the beam radius at the mirror. [m]
        zOff     - Initial guess of the z-offset. Generally unnecessary. [surfacemap.scaling]
        isCenter - 2D-list with booleans. isCenter[0] Determines if the center of the
                   sphere is to be fitted (True) or not (False, recommended). If the center is
                   fitted, then isCenter[1] determines if the weights (in case w!=None) are
                   centered around the fitted center (True) or centered around the center of
                   the data-grid (False, highly recommended).
                   
        if isCenter[0] == False
        Returns: Rc, zOff
        Rc       - Radius of curvature of the sphere removed from the mirror map. [m]
        zOff     - z-offset of the sphere removed. [surfacemap.scaling]
        isCenter[0] == True
        Returns Rc, zOff, center
        center   - [x0, y0] giving the coordinates of the center of the fitted sphere.
        
        Based on the file 'FT_remove_curvature_from_mirror_map.m' by Charlotte Bond.
        '''
        
        # z-value at centre of the data-grid. Serves as initial guess for deviation
        # from z=0, if no other first guess is given.
        if zOff is None:
            zOff = self.data[int(round(self.center[1])), int(round(self.center[0]))]

        # If fitting center of the sphere, four variables are fitted. Initial guess
        # of deviation from notNan-data-grid-center: (x0,y0) = (0,0).
        if isCenter[0]:
            p = [Rc0, zOff, 0.0, 0.0]
        # Otherwise two.
        else:
            p = [Rc0, zOff]

        # Grids with X,Y and r2 values. X and Y crosses zero in the center
        # of the xy-plane.
        X,Y = np.meshgrid(self.x, self.y)
        r2 = X**2 + Y**2
        
        # Cost-function to minimize.
        def costFunc(p):
            # If the center of the mirror is fitted, four variables are fitted.
            if isCenter[0]:
                Rc = p[0]
                zOff = p[1]
                x0 = p[2]
                y0 = p[3]
            # Otherwise 2 variables.
            else:
                Rc = p[0]
                zOff = p[1]
                x0 = 0
                y0 = 0

            Z = self.createSurface(Rc,X,Y,zOff,x0,y0)

            if w is None:
                # Mean squared difference between map and the created sphere.
                res = math.sqrt( ((self.data[self.notNan] -
                                   Z[self.notNan])**2).sum() )/self.notNan.sum()
            else:
                if isCenter[0] and isCenter[1]:
                    # Weights centered around fitting spehere center. May give weird
                    # results if the mirror deviates much from a sphere.
                    weight = (2/(math.pi*w**2))*np.exp(-2*( (X-x0)**2 + (Y-y0)**2 )/w**2)
                else:
                    # Weights centered around the center of the mirror xy-plane.
                    weight = (2/(math.pi*w**2))*np.exp(-2*r2/w**2)
                # Weighted mean squared difference between map and the created sphere.
                res = math.sqrt( ( weight[self.notNan]*
                                  ( (self.data[self.notNan] - Z[self.notNan])**2 )
                                  ).sum() ) \
                                  /weight[self.notNan].sum()
            return res

        # Using the simplex Nelder-Mead method. This is the same or very
        # similar to the method used in 'FT_remove_curvature_from_mirror_map.m',
        # but there are probably better methods to use.
        opts = {'xtol': 1.0e-5, 'ftol': 1.0e-9, 'maxiter': 10000, 'maxfev': maxfev, 'disp': False}
        out = minimize(costFunc, p, method='Nelder-Mead', options=opts)
        if not out['success']:
            msg = '  Warning: ' + out['message'].split('.')[0] + ' (nfev={:d}).'.format(out['nfev'])
            print(msg)
            
        # Assigning values to the instance variables
        self.RcRemoved = out['x'][0]
        if self.zOffset is None:
            self.zOffset = 0
        self.zOffset = self.zOffset + out['x'][1]

        # Equivalent Zernike (n=2,m=0) amplitude.
        R = self.find_radius(unit='meters')
        A20 = Rc2znm(self.RcRemoved,R)/self.scaling
        self.zernikeRemoved = (0,2,A20)

        # If center was fitted, assign new values to instance variable center, and
        # subtract the fitted sphere from the mirror map.
        if isCenter[0]:
            x0 = out['x'][2]
            y0 = out['x'][3]
            # Converts the deviation into a new absolut center in data points.
            self.center = (self.center[0] + x0/self.step_size[0],
                           self.center[1] + y0/self.step_size[1])
            # Creating fitted sphere
            Z = self.createSurface(self.RcRemoved, X, Y, self.zOffset, x0, y0)
            # Subtracting sphere from map
            self.data[self.notNan] = self.data[self.notNan]-Z[self.notNan]
            return self.RcRemoved, self.zOffset, self.center, A20
        # Subtracting fitted sphere from mirror map.
        else:
            # Creating fitted sphere
            Z = self.createSurface(self.RcRemoved,X,Y,self.zOffset)
            # Subtracting sphere from map
            self.data[self.notNan] = self.data[self.notNan]-Z[self.notNan]
            return self.RcRemoved, self.zOffset, A20

    def remove_curvature(self, method='zernike', w=None, zOff=None,
                         isCenter=[False,False], zModes = 'all'):
        '''
        Removes curvature from mirror map by fitting a sphere to the mirror map, or
        by convolving with second order Zernike polynomials and then extracting
        these polynomials with the obtained amplitudes. 

        Inputs: method, w, zOff, isCenter, zModes
        method   - 'zernike' (default) convolves second order Zernike polynomials with the
                   map, and then removes the polynomial with the obtained amplitude from the
                   surface map. Which of the three second order modes that are fitted is
                   determined by zMods (see below).
                 - 'sphere' fits a sphere to the mirror map. The fitting is Gaussian
                   weighted if w (see below) is specifed. 
        w        - Gaussian weighting parameter. The distance from the center where the
                   weigths have decreased by a factor of exp(-2) compared to the center.
                   Should preferrably be equal to the beam radius at the mirror. [m]
        zOff     - Initial guess of the z-offset. Only used together with method='sphere'.
                   Generally unnecessary. [surfacemap.scaling]
        isCenter - 2D-list with booleans. isCenter[0] Determines if the center of the
                   sphere is to be fitted (True) or not (False, recommended). If the center is
                   fitted, then isCenter[1] determines if the weights (in case w!=None) are
                   centered around the fitted center (True) or centered around the center of
                   the data-grid (False, highly recommended).
        zModes   - 'defocus' uses only Zernike polynomial (n,m) = (2,0), which has a
                   parabolic shape.
                   'astigmatism' uses Zernike polynomials (n,m) = (2,-2) and (n,m) = (2,2).
                   'all' uses both the defocus and the astigmatic modes.

        if method == 'zernike':
        Returns: Rc, znm
        Rc       - Approximate spherical radius of curvature removed.
        znm      - Dictionary specifying which Zernike polynomials that have been removed
                   as well as their amplitudes.
                   
        if method == 'sphere' and isCenter[0] == False:
        Returns: Rc, zOff
        Rc       - Radius of curvature of the sphere removed from the mirror map. [m]
        zOff     - z-offset of the sphere removed. [surfacemap.scaling]
        
        if method == 'sphere' and isCenter[0] == True:
        Returns Rc, zOff, center
        center   - [x0, y0] giving the coordinates of the center of the fitted sphere. 
        '''
        
        
        if method == 'zernike' or method == 'Zernike':
            Rc, znm = self.rmZernikeCurvs(zModes)
            return Rc, znm
        elif method == 'sphere' or method == 'Sphere':
            smap = deepcopy(self)
            # Using Zernike polynomial (m,n) = (0,2) to estimate radius of curvature.
            Rc0 = smap.rmZernikeCurvs('defocus')[0]
            if isCenter[0]:
                Rc, zOff, center, A20 = self.rmSphericalSurf(Rc0, w, zOff, isCenter)
                return Rc, zOff, self.center, A20
            else:
                Rc, zOff, A20 = self.rmSphericalSurf(Rc0, w, zOff, isCenter)
                return Rc, zOff, A20

            
    def recenter(self):
        '''
        Setting center of mirror to be in the middle of the notNan field. 
        
        Based on 'FT_recenter_mirror_map.m' by Charlotte Bond.
        '''
        # Row and column indices with non-NaN elements
        rIndx, cIndx = self.notNan.nonzero()
        # Finding centres
        x0 = float(cIndx.sum())/len(cIndx)
        y0 = float(rIndx.sum())/len(rIndx)
        self.center = tuple([x0,y0])
        self.xyOffset = (.0,.0)
        return self.center
        # -------------------------------------------------
        
    def crop(self, radius=None):
        '''
        If no radius is given it finds the NaN elements closest to the map
        center, and sets all elements further out to NaN. Then cropping
        and recentering. If radius is given, it sets all elements outside
        the radius value to NaN, cropping, and recentering.

        Inputs: radius
        radius - crops removes elements outside this value [m].
        
        Based on Charlotte Bond's simtools function
        'FT_remove_elements_outside_map.m'.
        '''
        # Saves the xyOffset to restore it afterwards
        xyOffset = self.xyOffset
        self.recenter()
        if radius is None:
            # Arrays of row and column indices where data is NaN.
            r,c = np.where(self.notNan == False)
            # Sets the map center to index [0,0]
            x = c-self.center[0]
            y = r-self.center[1]
            # Minimum radius squared measured in data points.
            r2 = (x**2+y**2).min()
        else:
            # Requires the same stepsize in x and y.
            r2 = (radius/self.step_size[0])**2
            
        # Grid with the same dimensions as the map.
        X,Y = np.meshgrid(self.x, self.y)
        # Grid containing radial distances squared measured from the center.
        R2 = (X/self.step_size[0])**2 + (Y/self.step_size[1])**2
        # Matrix with True=='outside mirror surface'
        outs = R2>=r2
        # Setting notNan-matrix to false outside the mirror.
        self.notNan[outs] = False
        # Setting data matrix to 0 outside the mirror.
        self.data[outs] = 0.0

        # Removing unnecessary data points.
        # --------------------------------------
        # Radius inside data is kept. Don't know why R is set this way,
        # but trusting Dr. Bond.
        R = round(math.ceil(2.0*math.sqrt(r2)+6.0)/2.0)
        if (2*R<self.data.shape[0] and 2*R<self.data.shape[1] and
            R<self.center[0] and R<self.center[1] and
            R<(self.size[0]-self.center[0]) and R<(self.size[1]-self.center[1])):
            
            x0 = round(self.center[0])
            y0 = round(self.center[1])

            self.data = self.data[int(y0-R):int(y0+R+1), int(x0-R):int(x0+R+1)]
            self.notNan = self.notNan[int(y0-R):int(y0+R+1), int(x0-R):int(x0+R+1)]

        # Centering the new cropped map and restoring offset
        self.recenter()
        self.xyOffset = xyOffset
        
    
    def createSurface(self,Rc,X,Y,zOffset=0,x0=0,y0=0,xTilt=0,yTilt=0,isPlot=False):
        '''
        Creating surface.
        
        Inputs: Rc, X, Y, zOffset=0, x0=0, y0=0, xTilt=0, yTilt=0, isPlot=False
        Rc      - Radius of curvature, and center of sphere on z-axis in case zOffset=0. [m]
        X, Y    - Matrices of x,y-values generated by 'X,Y = numpy.meshgrid(xVec,yVec)'. [m]
        zOffset - Surface center offset from 0. Positive value gives the surface center
                  positive z-coordinate. [self.scaling]
        x0,y0   - The center of the sphere in the xy-plane.
        x/yTilt - Tilt of x/y-axis [rad]. E.g. xTilt rotates around the y-axis.
        isPlot  - Plot or not [boolean]. Not recommended when used within optimization
                  algorithm.

        Returns: Z
        Z       - Matrix of surface height values. [self.scaling]
        
        Based on Simtools function 'FT_create_sphere_for_map.m' by Charlotte Bond.
        '''
        
        if Rc == np.inf:
            Rc = 0
        
        # Adjusting for tilts and offset
        Z = zOffset + (X*np.tan(xTilt) + Y*np.tan(yTilt))/self.scaling
        
        # Adjusting for spherical shape.
        if Rc != 0 and Rc is not None:
            Z = Z + (Rc - np.sign(Rc)*np.sqrt(Rc**2 - (X-x0)**2 - (Y-y0)**2))/self.scaling
        
        if isPlot:
            import pylab
            pylab.clf()
            fig = pylab.figure()
            axes = pylab.pcolormesh(X, Y, Z) #, vmin=zmin, vmax=zmax)
            cbar = fig.colorbar(axes)
            cbar.set_clim(Z.min(), Z.max())
            pylab.title('Z_min = ' + str(Z.min()) + '   Z_max = ' + str(Z.max()))
            pylab.show()
        return Z
    
    def createPolarGrid(self):
        '''
        Creating polar grid from the rectangular x and y coordinates in
        surfacemap. 

        Returns: rho, phi
        rho  - matrix with radial distances from centre of mirror map.
        phi  - matrix with polar angles.
        '''
        X,Y = np.meshgrid(self.x,self.y)
        phi = np.arctan2(Y,X)
        rho = np.sqrt(X**2 + Y**2)
        return rho, phi

    def preparePhaseMap(self, w=None, xyOffset=None, verbose=False):
        '''
        Prepares the phase mirror map for being used in Finesse. The map is cropped, and
        the curvature, offset, and tilts are removed, in that order.

        Input: w, xyOffset, verbose
        w        - radius of Gaussian weights [m]. If specified, surfaces are fitted
                   to the mirror map when removing curvature and tilt. Otherise, the
                   whole mirror map is treated the same, thus covolutions with Zernike
                   polynomials are used instead.
        xyOffset - Beam center offset compared to the mirror center.
        verbose  - True makes a chatty method, while False only prints the end results.

        Output: amap
        amap     - Aperture map (absorption map), i.e. zeros within the map surface and
                   ones outside it.
                   
        Based on Charlotte Bond's Simtools function 'FT_prepare_phase_map_for_finesse.m'.
        '''
        if verbose:
            print('--------------------------------------------------')
            print('Preparing phase map for Finesse...')
            print('--------------------------------------------------')
        w_max = self.find_radius(method='min', unit='meters')
        if w is None:
            if verbose:
                print(' No weighting used')
        elif w >= w_max:
            if verbose:
                print(' Weighting radius ({:.3f} cm) larger than mirror radius ({:.3f} cm).'.format(w*100, w_max*100))
                print(' Thus, no weighting used.')
            w = None
        elif verbose:
            print(' Gaussian weights used with radius: {:.2f} cm'.format(w*100.0))
        if verbose:
            print(' (rms, avg) = ({:.3e}, {:.3e}) nm'.format(self.rms(w),self.avg(w)))
            print('--------------------------------------------------')
        
            print(' Centering...')
        self.recenter()
        if verbose:
            print('  New center (x0, y0) = ({:.2f}, {:.2f})'.
                  format(self.center[0],self.center[1]))
        if verbose:
            print(' Cropping map...')
        before = self.size
        self.crop()
        if verbose:
            print('  Size (rows, cols): ({:d}, {:d}) ---> ({:d}, {:d})'.
                  format(before[1], before[0], self.size[1], self.size[0]))
            print('  New center (x0, y0) = ({:.2f}, {:.2f})'.
                  format(self.center[0],self.center[1]))
            print('  (rms, avg) = ({:.3e}, {:.3e}) nm'.
                  format(self.rms(w),self.avg(w)))

        # Radius of mirror in xy-plane.
        R = self.find_radius(unit='meters')
        
        if verbose:
            print(' Removing curvatures...')
        # --------------------------------------------------------
        if w is None:
            Rc, znm = self.remove_curvature(method='zernike', zModes = 'defocus')
            if verbose:
                print('  Removed Z20 with amplitude A20 = {:.2f} nm'.format(znm['02'][2]))
                print('  Equivalent Rc = {:.2f} m'.format(Rc))
        else:
            Rc, zOff, A20 = self.remove_curvature(method='sphere', w=w)
            # Adding A20 to zOff (calling it A0 now) because: Z(n=2,m=0) = A20*(r**2 - 1)
            A0 = zOff+A20
            if verbose:
                print('  Removed Rc = {0:.2f} m'.format(Rc) )
                print('  Equivalent Z(n=2,m=0) amplitude A20 = {0:.2f} nm'.format(A20))

        if verbose:
            print('  (rms, avg) = ({:.3e}, {:.3e}) nm'.format(self.rms(w),self.avg(w)))
            
            print(' Removing offset...')
        # --------------------------------------------------------

        zOff = self.removeOffset(w)
        if w is None:
            if verbose:
                print('  Removed Z00 with amplitude A00 = {:.2f} nm'.format(zOff) )
        else:
            A0 = A0 + zOff
            if verbose:
                print('  Removed z-offset (A00) = {0:.3f} nm'.format(zOff))
        if verbose:
            print('  (rms, avg) = ({:.3e}, {:.3e}) nm'.format(self.rms(w),self.avg(w)))
        
            print(' Removing tilts...')
        # --------------------------------------------------------
        if w is None:
            A1,xbeta,ybeta = self.rmTilt(method='zernike')
            if verbose:
                for k in range(len(A1)):
                    print('  Removed Z1{:d} with amplitude A1{:d} = {:.2f} nm'.
                          format(2*k-1,2*k-1,A1[k]) )
                print('  Equivalent tilts in radians: ')
                print('   xbeta = {:.2e} rad'.format(xbeta))
                print('   ybeta = {:.2e} rad'.format(ybeta))
        else:
            A1,xbeta,ybeta,zOff = self.rmTilt(method='fitSurf', w=w)
            A0 = A0 + zOff
            if verbose:
                print('  Tilted surface removed:')
                print('   xbeta    = {:.2e} rad'.format(xbeta))
                print('   ybeta    = {:.2e} rad'.format(ybeta))
                print('   z-offset = {:.2e} nm'.format(zOff))
                print('  Equivalent Zernike amplitudes:')
                print('   A(1,-1) = {:.2f} nm'.format(A1[0]))
                print('   A(1, 1) = {:.2f} nm'.format(A1[1]))
        if verbose:
            print('  (rms, avg) = ({:.3e}, {:.3e}) nm'.format(self.rms(w),self.avg(w)))

        if w is not None:
            # Adding the total offset removed to zernikeRemoved (just for the record)
            self.zernikeRemoved = (0,0,A0)
            if verbose:
                print(' Equivalent Z00 amplitude from accumulated z-offsets, A00 = {:.2f}'.format(A0))
                
        if verbose:
            print(' Creating aperture map...')
            # --------------------------------------------------------
        amap = aperturemap(self.name, self.size, self.step_size,
                           self.find_radius(method='min',unit='meters'), self.center)
        if verbose:
            print('  Aperture map created.')
                                       
        if xyOffset is not None:
            if verbose:
                print(' Offsetting mirror center in the xy-plane...')
            # --------------------------------------------------------
            self.center = (self.center[0] + xyOffset[0]/self.step_size[0],
                           self.center[1] + xyOffset[1]/self.step_size[1])
            self.xyOffset = xyOffset
            if verbose:
                print('  New mirror center (x0, y0) = ({:.2f}, {:.2f})'.format(self.center[0],self.center[1]))
        else:
            self.xyOffset = (0,0)
            
        if verbose:
            print(' Writing maps to file...')
        # --------------------------------------------------------
        filename = self.name + '_finesse.txt'
        self.write_map(filename)
        if verbose:
            print('  Phase map written to file {:s}'.format(filename))
        filename = amap.name + '_aperture.txt'
        amap.write_map(filename)
        if verbose:
            print('  Aperture map written to file {:s}'.format(filename))
        if verbose:
            print(' Writing result information to file...')
            # --------------------------------------------------------
        filename = self.name + '_finesse_info.txt'
        self.writeMapPrepResults(filename,w)
        if verbose:
            print('  Result written to file {:s}'.format(filename))
            
            
        # Printing results to terminal
        print('--------------------------------------------------')
        print('Phase map prepared!' )
        print('Name: {:s}'.format(self.name))
        print('Type: {:s}'.format(self.type))
        print('--------------------------------------------------')

        if w is None:
            print('Weights:   None')
        else:
            print('Weights:   r     = {:.2f} cm'.format(w*100))
        print('Radius:    R     = {:.2f} cm'.format(self.find_radius(unit='meters')*100))
        print('Offset:    A00   = {:.2f} nm'.format(self.zernikeRemoved['00'][2]))
        print('Tilt:      A1-1  = {:.2f} nm'.format(self.zernikeRemoved['-11'][2]))
        print('           A11   = {:.2f} nm,  or'.format(self.zernikeRemoved['11'][2]))
        print('           xbeta = {:.2e} rad'.format(self.betaRemoved[0]))
        print('           ybeta = {:.2e} rad'.format(self.betaRemoved[1]))
        print('Curvature: A20   = {:.2f} nm,  or'.format(self.zernikeRemoved['02'][2]))
        print('           Rc    = {:.2f} m'.format(self.RcRemoved))
        print('xy-offset: x0    = {:.2f} mm'.format(self.xyOffset[0]*1000))
        print('           y0    = {:.2f} mm'.format(self.xyOffset[1]*1000))
        print('Stats:     rms   = {:.3e} nm'.format(self.rms(w)))
        print('           avg   = {:.3e} nm'.format(self.avg(w)))
        print('--------------------------------------------------')
        print()

        return amap
    
    def writeMapPrepResults(self, filename, w):
        '''
        Writing result information to file. The same info that is written to the
        terminal.
        '''
        import time
        with open(filename,'w') as mapfile:
            mapfile.write('---------------------------------------\n')
            mapfile.write('Map:  {:s}\n'.format(self.name))
            mapfile.write('Date:  {:s}\n'.format(time.strftime("%d/%m/%Y  %H:%M:%S")))
            mapfile.write('---------------------------------------\n')
            if w is None:
                mapfile.write('Weights:   None\n')
            else:
                mapfile.write('Weights:   r     = {:.2f} cm\n'.format(w*100))
            mapfile.write('Diameter:  D     = {:.2f} cm\n'.format(self.find_radius(unit='meters')*200.0))
            mapfile.write('Offset:    A00   = {:.2f} nm\n'.format(self.zernikeRemoved['00'][2]))
            mapfile.write('Tilt:      A1-1  = {:.2f} nm\n'.format(self.zernikeRemoved['-11'][2]))
            mapfile.write('           A11   = {:.2f} nm,  or\n'.format(self.zernikeRemoved['11'][2]))
            mapfile.write('           xbeta = {:.2e} rad\n'.format(self.betaRemoved[0]))
            mapfile.write('           ybeta = {:.2e} rad\n'.format(self.betaRemoved[1]))
            mapfile.write('Curvature: A20   = {:.2f} nm,  or\n'.format(self.zernikeRemoved['02'][2]))
            mapfile.write('           Rc    = {:.2f} m\n'.format(self.RcRemoved))
            mapfile.write('xy-offset: x0    = {:.2f} mm\n'.format(self.xyOffset[0]*1000))
            mapfile.write('           y0    = {:.2f} mm\n'.format(self.xyOffset[1]*1000))
            
            # mapfile.write('Offset (Z00):  {:.2f}'.format(self.zernikeRemoved['00'][2]))
          
        
    def removeOffset(self, r=None):
        '''
        Based on Simtools function 'FT_remove_offset_from_mirror_map.m' by Charlotte Bond.
        '''
        if r is None:
            offset = self.rmZernike(0,0)
        else:
            rho = self.createPolarGrid()[0]
            inside = rho<r
            inside_notNan = np.zeros(inside.shape,dtype=bool)
            inside_notNan[inside] = self.notNan[inside]

            offset = self.data[inside_notNan].sum()/inside_notNan.sum()
            self.data[self.notNan] = self.data[self.notNan]-offset
        
        return offset

    def rmZernike(self, n, m='all'):
        '''
        Removes Zernike polynomials from mirror map.

        Input: n, m
        n  - Radial degree of Zernike polynomial.
        m  - Azimuthal degree of Zernike polynomial. If set to 'all', all polynomials
             of radial degree n is removed. Can also be an integer, or a list
             specifying the Azimuthal degrees that will be removed. Must be consistent
             with the chosen radial degree.

        Returns: A
        A  - List of amplitudes of the removed Zernike polynomials.
        '''
        if m == 'all':
            m = list(range(-n,n+1,2))

        R = self.find_radius(unit='meters')
        A = self.zernikeConvol(n,m)
        rho, phi= self.createPolarGrid()
        rho = rho/R

        if isinstance(m, list):
            for k in range(len(m)):
                Z = zernike(m[k], n, rho, phi)
                self.data[self.notNan] = self.data[self.notNan]-A[k]*Z[self.notNan]
                self.zernikeRemoved = (m[k], n, A[k])
        else:
            Z = zernike(m, n, rho, phi)
            self.data[self.notNan] = self.data[self.notNan]-A*Z[self.notNan]
            self.zernikeRemoved = (m, n, A)
            
        return A
        
class mergedmap:
    """
    A merged map combines multiple surfaces map to form one. Such a map can be used
    for computations of coupling coefficients but it cannot be written to a file to 
    be used with Finesse. For this you must output each map separately.
    
    """
    
    def __init__(self, name, size, center, step_size, scaling):
        
        self.name = name
        self.center = center
        self.step_size = step_size
        self.scaling = scaling
        self.__interp = None
        self._rom_weights = None
        self.__maps = []
        self.weighting = None
        
    def addMap(self, m):
        self.__maps.append(m)
    
    @property
    def center(self):
        return self.__center
    
    @center.setter
    def center(self, value):
        self.__center = value
        self.__interp = None
    
    @property
    def type(self):
        hasR = False
        hasT = False
        
        _type = ""
        
        for m in self.__maps:
            if "reflection" in m.type: hasR = True
            
            if "transmission" in m.type: hasT = True
            
            if "both" in m.type:
                hasR = True
                hasT = True
        
        if hasR and not hasT: _type += "reflection "
        elif hasR and not hasT: _type += "transmission "
        elif hasR and hasT: _type += "both "
        
        return _type
        
    @property
    def step_size(self):
        return self.__step_size
    
    @step_size.setter
    def step_size(self, value):
        self.__step_size = value
        self.__interp = None

    @property
    def scaling(self):
        return self.__scaling
    
    @scaling.setter
    def scaling(self, value):
        self.__scaling = value
        self.__interp = None
    
    @property
    def x(self):
        return self.step_size[0] * (np.array(range(0, self.size[0])) - self.center[0])
        
    @property
    def y(self):
        return self.step_size[1] * (np.array(range(0, self.size[1])) - self.center[1])

    @property
    def size(self):
        return self.__maps[0].data.shape
            
    @property
    def offset(self):
        return np.array(self.step_size)*(np.array(self.center) - (np.array(self.size)-1)/2.0)
    
    @property
    def ROMWeights(self):
        return self._rom_weights
    
    def z_xy(self, wavelength=1064e-9, direction="reflection_front", nr1=1.0, nr2=1.0):
        
        z_xy = np.ones(self.size, dtype=np.complex128)
        
        for m in self.__maps:
            z_xy *= m.z_xy(wavelength=wavelength, direction=direction, nr1=nr1, nr2=nr2)
            
        if self.weighting is None:
            return z_xy
        else:
            return z_xy * self.weighting
        
    def generateROMWeights(self, EIxFilename, EIyFilename=None, verbose=False, interpolate=False, newtonCotesOrder=8, nr1=1, nr2=1):
        if interpolate == True:
            # Use EI nodes to interpolate if we
            with open(EIxFilename, 'rb') as f:
                EIx = pickle.load(f)

            if EIyFilename is None:
                EIy = EIx
            else:
                with open(EIyFilename, 'rb') as f:
                    EIy = pickle.load(f)

            x = EIx.x
            x.sort()
            nx = np.unique(np.hstack((x, -x[::-1])))
        
            y = EIy.x
            y.sort()
            ny = np.unique(np.hstack((y, -y[::-1])))
            
            self.interpolate(nx, ny)
        
        w_refl_front, w_refl_back, w_tran_front, w_tran_back = (None, None, None, None)
        
        if "reflection" in self.type or "both" in self.type:
            w_refl_front = makeWeightsNew(self, EIxFilename, EIyFilename,
                                      verbose=verbose, newtonCotesOrderMapWeight=newtonCotesOrder,
                                      direction="reflection_front")
            
            w_refl_front.nr1 = nr1
            w_refl_front.nr2 = nr2
            
            w_refl_back = makeWeightsNew(self, EIxFilename, EIyFilename,
                                      verbose=verbose, newtonCotesOrderMapWeight=newtonCotesOrder,
                                      direction="reflection_back")
            
            w_refl_back.nr1 = nr1
            w_refl_back.nr2 = nr2

        if "transmission" in self.type or "both" in self.type:                                      
            w_tran_front = makeWeightsNew(self, EIxFilename, EIyFilename,
                                      verbose=verbose, newtonCotesOrderMapWeight=newtonCotesOrder,
                                      direction="transmission_front")

            w_refl_front.nr1 = nr1
            w_refl_front.nr2 = nr2
                                            
            w_tran_back  = makeWeightsNew(self, EIxFilename, EIyFilename,
                                      verbose=verbose, newtonCotesOrderMapWeight=newtonCotesOrder,
                                      direction="transmission_back")
            
            w_refl_back.nr1 = nr1
            w_refl_back.nr2 = nr2
            
        self._rom_weights = MirrorROQWeights(w_refl_front, w_refl_back, w_tran_front, w_tran_back)
        
        return self._rom_weights

    def interpolate(self, nx, ny, **kwargs):
        """
        Interpolates all the maps that are used to fc
        
        Uses scipy.interpolate.interp2d and any keywords arguments are
        passed on to it, thus settings like interpolation type and
        fill values can be set.
        
        The range of nx and ny must contain the value zero so that the
        center point of the map can be set.
        """

        for m in self.__maps:
            m.interpolate(nx, ny)

    def plot(self, mode="absorption", show=True, clabel=None, xlim=None, ylim=None, wavelength=1064e-9, isBlock=False):
        
        import matplotlib
        import matplotlib.pyplot as plt
        
        if xlim is not None:
            _x = np.logical_and(self.x<=max(xlim)/100.0, self.x>=min(xlim)/100.0)
            xmin = np.min(np.where(_x == True))
            xmax = np.max(np.where(_x == True))
        else:
            xmin = 0
            xmax = len(self.x)-1
            xlim = [self.x.min()*100, self.x.max()*100]
    
        if ylim is not None:
            _y = np.logical_and(self.y<=max(ylim)/100.0, self.y>=min(ylim)/100.0)
            ymin = np.min(np.where(_y == True))
            ymax = np.max(np.where(_y == True))
        else:
            ymin = 0
            ymax = len(self.y)-1
            ylim = [self.y.min()*100, self.y.max()*100]

        if mode == "absorption":
            # plots how much of field is absorbed
            data = 1-np.abs(self.z_xy())
        elif mode == "meter":
            # plot the phase in terms of meters of displacement
            k = 2*np.pi/wavelength
            data = np.angle(self.z_xy()) / (2*k)
            
        zmin = data[xmin:xmax,ymin:ymax].min()
        zmax = data[xmin:xmax,ymin:ymax].max()

        # 100 factor for scaling to cm
        xrange = 100*self.x
        yrange = 100*self.y

        fig = plt.figure()
        axes = plt.pcolormesh(xrange, yrange, data, vmin=zmin, vmax=zmax)
        plt.xlabel('x [cm]')
        plt.ylabel('y [cm]')

        if xlim is not None: plt.xlim(xlim)
        if ylim is not None: plt.ylim(ylim)

        plt.title('Merged map {0}, mode {1}'.format(self.name, mode))

        cbar = fig.colorbar(axes)
        cbar.set_clim(zmin, zmax)
        
        if clabel is not None:
            cbar.set_label(clabel)
    
        if show:
            plt.show(block=isBlock)
        
        return fig

class aperturemap(surfacemap):
    
    def __init__(self, name, size, step_size, R, center=None):
        if center is None:
            center = (np.array(size)+1)/2.0
        surfacemap.__init__(self, name, "absorption both", size, center, step_size, 1)
        self.R = R
        
    @property
    def R(self):
        return self.__R
    
    @R.setter
    def R(self, value):
        self.__R = value
    
        xx, yy = np.meshgrid(self.x, self.y)
        
        radius = np.sqrt(xx**2 + yy**2)
        
        self.data = np.zeros(self.size[::-1])
        self.data[radius > self.R] = 1.0
        
        
class curvedmap(surfacemap):
    """
    Example:
        import pykat
        from pykat.optics.knm import *
        from pykat.optics.maps import *

        N = 501
        dx = 1/float(N)
        Roc = 10e3
    
        m = curvedmap("test", (N,N), (dx,dx), RoC)
    """
    def __init__(self, name, size, step_size, Rc):
        surfacemap.__init__(self, name, "phase reflection", size, (np.array(size)+1)/2.0, step_size, 1e-6)
        
        self.Rc = Rc
            
    @property
    def Rc(self):
        return self.__Rc
    
    @Rc.setter
    def Rc(self, value):
        self.__Rc = float(value)
    
        xx, yy = np.meshgrid(self.x, self.y)
        
        Rsq = xx**2 + yy**2

        if np.any(self.Rc**2-Rsq[self.notNan] < 0):
            raise BasePyKatException("Invalid curvature Rc, must be bigger than radius of the mirror")
        else:
            self.data[self.notNan] = (self.Rc - math.copysign(1.0, self.Rc) * np.sqrt(self.Rc**2 - Rsq[self.notNan]))/ self.scaling

class tiltmap(surfacemap):
    """
    To create a tiltmap, plot it and write it to a file to use with Finesse:
        
        tilts = (1e-6, 1e-8) # tilt in (x, y) radians\
        dx = 1e-4
        L = 0.2
        N = L/dx
        
        tmap = tiltmap("tilt", (N, N), (dx,dx), tilts)
        tmap.plot()
        tmap.write_map("mytilt.map")
    """
    
    def __init__(self, name, size, step_size, tilt):
        surfacemap.__init__(self, name, "phase reflection", size, (np.array(size)+1)/2.0, step_size, 1e-9)
        self.tilt = tilt
        
    @property
    def tilt(self):
        return self.__tilt
    
    @tilt.setter
    def tilt(self, value):
        self.__tilt = value
        
        xx, yy = np.meshgrid(self.x, self.y)
        
        self.data = (yy * self.tilt[1] + xx * self.tilt[0])/self.scaling
        

class zernikemap(surfacemap):
	def __init__(self, name, size, step_size, radius, scaling=1e-9):
		surfacemap.__init__(self, name, "phase reflection", size, (np.array(size)+1)/2.0, step_size, scaling)
		self.__zernikes = {}
		self.radius = radius
		
	@property
	def radius(self): return self.__radius

	@radius.setter
	def radius(self, value, update=True):
		self.__radius = float(value)
		if update: self.update_data()

	def setZernike(self, m, n, amplitude, update=True):
		self.__zernikes["%i%i" % (m, n)] = (m,n,amplitude)
		if update: self.update_data()

	def update_data(self):
		X,Y = np.meshgrid(self.x, self.y)
		R = np.sqrt(X**2 + Y**2)
		PHI = np.arctan2(Y, X)

		data = np.zeros(np.shape(R))

		for i in self.__zernikes.items():
			data += i[1][2] * zernike(i[1][0], i[1][1], R/self.radius, PHI)

		self.data = data
	
			

  
def read_map(filename, mapFormat='finesse', scaling=1.0e-9, mapType='phase', field='both'):
    '''
    Reads surface map files and returns a surfacemap object.
    
    filename  - name of surface map file.
    mapFormat - 'finesse', 'ligo', 'zygo', 'metroPro' (binary).
    scaling   - scaling of surface height of the mirror map [m].
    '''
    
    # Reads finesse mirror maps.
    if mapFormat.lower() == 'finesse':
        
        g = lambda x: float(x)
        with open(filename, 'r') as f:
        
            f.readline()
            name = f.readline().split(':')[1].strip()
            maptype = f.readline().split(':')[1].strip()
            size = tuple(map(g, f.readline().split(':')[1].strip().split()))
            center = tuple(map(g, f.readline().split(':')[1].strip().split()))
            step = tuple(map(g, f.readline().split(':')[1].strip().split()))
            scaling = float(f.readline().split(':')[1].strip())
        
        data = np.loadtxt(filename, dtype=np.float64,ndmin=2,comments='%')
         
        return surfacemap(name, maptype, size, center, step, scaling, data)
        
    elif mapFormat.lower() == 'ligo' or mapFormat.lower() == 'zygo':
        '''
        Reads zygo and ligo mirror maps. Based on Charlotte Bond's Simtools
        functions 'FT_read_zygo_map.m' and 'FT_read_ligo_map.m'.
        '''
        maptype = ' '.join([mapType, field])
        if mapFormat == 'ligo':
            isLigo = True
            # Unused, but need to set it to call the function.
            isAscii = False
            # Remove '_asc.dat' for output name
            name = filename.split('_')
            name = '_'.join(name[:-1])
        else:
            isLigo = False
            tmp = filename.split('.')
            fileFormat = tmp[-1].strip()
            name = '.'.join(tmp[:-1])
            if fileFormat == 'asc':
                isAscii = True
            else:
                isAscii = False
                
        # Reading the zygo/ligo surface map file. Note that the intensity data
        # iData is available here too when ascii zygo files are read, but at
        # the moment it's not used. 
        out = readZygoLigoMaps(filename, isLigo, isAscii)
        if len(out) == 3:
            hData, data, iData = out
        else:
            hData, data = out
        
        # Matrix with True where data element is NaN.
        isNan = np.isnan(data)
        # Matrix with True where data element is not NaN.
        notNan = isNan==False
        # Setting NaN-values to zero
        data[isNan] = 0
        # Scaling to the scaling chosesen as input.
        data[notNan] = (hData['hScale']/scaling)*data[notNan]
        smap = surfacemap(name, maptype, hData['size'], hData['center'],
                          hData['stepSize'], scaling, data, notNan)
        smap.recenter()
        
        return smap

        
    if mapFormat == 'metroPro':
        # Remove directory part and '.dat' for name of map
        name = filename.split('/')[-1].split('.')[0]
        # Reading file
        data, hData = readMetroProData(filename)
        # xy-stepsize in meters
        stepSize = (hData['cameraRes'], hData['cameraRes'])
        # Temporary center
        center = ((len(data[0,:])-1)/2.0,(len(data[:,0])-1)/2.0)
        # Keeping track of nan's
        isNan = np.isnan(data)
        notNan = isNan==False
        data[isNan] = 0
        # Type of map
        maptype = ' '.join([mapType, field])
        # Unnecessary size, remove from class later
        size = data.shape[::-1]
        # Creating surfacemap-instance
        smap = surfacemap(name, maptype, size, center, stepSize, 1.0, data, notNan)
        smap.rescaleData(scaling)
        smap.recenter()
        return smap


def readZygoLigoMaps(filename, isLigo=False, isAscii=True):
    '''
    Converts raw zygo and ligo mirror maps to the finesse
    format. Based on the matlab scripts 'FT_read_zygo_map.m' and
    'FT_read_ligo_map.m'.
    '''
    g = lambda x: float(x)
    hData = {}
    with open(filename, 'r') as f:
        # Reading header of LIGO-map
        # ------------------------------------------------------
        # Skip first two lines
        for k in range(2):
            f.readline()
        # If zygo-file, and ascii format, there is intensity
        # data. Though, the Ligo files I have seen are also
        # zygo-files, so maybe we should extract this data
        # from these too?
        line = f.readline()
        if not isLigo and isAscii:
            iCols = float(line.split()[2])
            iRows = float(line.split()[3])
            
        line = f.readline().split()
        
        if isLigo:
            y0 = float(line[0])
            x0 = float(line[1])
            rows = float(line[2])
            cols = float(line[3])
        else:
            y0 = float(line[1])
            x0 = float(line[0])
            rows = float(line[3])
            cols = float(line[2])

        # Skipping three lines
        for k in range(3):
            f.readline()
        line = f.readline().split()

        # Scaling factors
        # ----------------------------------------------
        # Interfeometric scaling factor
        S = float(line[1])
        # Wavelength
        lam = float(line[2])
        # Obliquity factor
        O = float(line[4])
        # ----------------------------------------------
        # Physical step size in metres
        if line[6] != 0:
            xstep = float(line[6])
            ystep = float(line[6])
        else:
            xstep = 1.0
            ystep = 1.0

        # Skipping two lines
        for k in range(2):
            f.readline()
        line = f.readline().split()

        # Resolution of phase data points, 1 or 0.
        phaseRes = float(line[0])
        if phaseRes == 0:
            R = 4096
        elif phaseRes == 1:
            R = 32768

        if not isLigo and not isAscii:
            # zygo .xyz files give phase data in microns.
            hData['hScale'] = 1.0e-6
            
        else:
            # zygo .asc and ligo-files give phase data in
            # internal units. To convert to m use hScale
            # factor.
            hData['hScale'] = S*O*lam/R

        # Skipping three lines
        for k in range(3):
            f.readline()
            
        if not isLigo:
            if not isAscii:
                # For the zygo .xyz-format.
                k = 0
                data = np.zeros(int(rows*cols))
                totRuns = cols*rows
            
                # Read data
                while k < cols*rows:
                    line = f.readline().split()
                    # Setting 'no data' to NaN
                    if len(line)==4:
                        data[k]=np.nan
                    # Reading value
                    elif len(line)==3:
                        data[k] = float(line[2])
                    # Otherwise stop reading
                    else:
                        k = cols*rows
                    k+=1
            else:
                # Skipping one line
                line = f.readline()
                # Reading intensity data
                iData = np.array([])
                line = f.readline().split()
                while line[0] != '#':
                    iData = np.append(iData, [x for x in map(g,line)])
                    line = f.readline().split()
                # Reshaping intensity data
                iData = iData.reshape(int(iRows), int(iCols)).transpose()
                iData = np.rot90(iData)
        else:
            # Skipping one line
            f.readline()
            # Skipping lines until '#' is found.
            while f.readline()[0] != '#':
                pass

        if isLigo or (not isLigo and isAscii):
            # Reading phase data
            # ----------------------------------------------
            # Array with the data
            data = np.array([])
            # Reading data until next '#' is reached.
            line = f.readline().split()
            while line[0] != '#':
                data = np.append(data, [x for x in map(g,line)])
                line = f.readline().split()
            # ----------------------------------------------


    if isLigo:
        # Setting all the points outside of the mirror
        # surface to NaN. These are given a large number
        # in the file. 
        data[data == data[0]] = np.nan

        # Reshaping into rows and columns
        data = data.reshape(int(cols),int(rows)).transpose()
        # Pretty sure that the lines below can be done in
        # less operations, but it's quick as it is.
        # ----------------------------------------------
        # Flipping right and left
        data = np.fliplr(data)
        # Rotating 90 degrees clockwise 
        data = np.rot90(data,-1)
        # Flipping right and left
        data = np.fliplr(data)
        # ----------------------------------------------
    else:
        if isAscii:
            # Setting all the points outside of the mirror
            # surface to NaN. These are given a large number
            # in the file. 
            data[data >= 2147483640] = np.nan
        # Reshaping into rows and columns.
        data = data.reshape(int(rows),int(cols)).transpose()
        # Rotating to make (0,0) be in bottom left
        # corner. 
        data = np.rot90(data)

    hData['size'] = data.shape[::-1]

    # Center according to file
    hData['center'] = tuple([x0,y0])
    # xy-stepsize
    hData['stepSize'] = tuple([xstep,ystep])

    if not isLigo and isAscii:
        return hData, data, iData
    else:
        return hData, data



def readMetroProData(filename):
    '''
    Reading the metroPro binary data files.

    Translated from Hiro Yamamoto's 'LoadMetroProData.m'.
    '''
    f = BinaryReader(filename)
    # Read header
    hData = readHeaderMP(f)
    if hData['format'] < 0:
        print('Error: Format unknown to readMetroProData()\nfilename: {:s}'.format(filename))
        return 0
    # Read phase map data
    # Skipping header and intensity data
    f.seek(hData['size']+hData['intNBytes'])
    # Reading data
    dat = f.read('int32',size=hData['Nx']*hData['Ny'])
    # Marking unmeasured data as NaN
    dat[dat >= hData['invalid']] = np.nan
    # Scale data to meters
    dat = dat*hData['convFactor']
    # Reshaping into Nx * Ny matrix
    dat = dat.reshape(hData['Ny'], hData['Nx'])
    # Flipping up/down, i.e., change direction of y-axis.
    dat = dat[::-1,:]
    # Auxiliary data to return
    dxy = hData['cameraRes']

    # Unnecessary
    #x1 = dxy * np.arange( -(len(dat[0,:])-1)/2, (len(dat[0,:])-1)/2 + 1)
    #y1 = dxy * np.arange( -(len(dat[:,0])-1)/2, (len(dat[:,1])-1)/2 + 1)

    return dat, hData #, x1, y1


def readHeaderMP(f):
    '''
    Reads header of the metroPro format files.
    
    Translated from the readerHeader() function within the
    'LoadMetroProData.m' function written by Hiro Yamamoto.
    '''
    hData = {}
    hData['magicNum'] = f.read('uint32')
    hData['format'] = f.read('int16')
    hData['size'] = f.read('int32')
    # Check if the magic string and format are known ones.
    if not (hData['format'] >= 1 and hData['format']<=3 and
            hData['magicNum']-hData['format'] == int('881B036E',16)):
        hData['format'] = -1
    # Read necessary data
    hData['invalid'] = int('7FFFFFF8',16)
    f.seek(60)
    # Intensity data, which we will skip over.
    hData['intNBytes'] = f.read('int32')
    # Top-left coordinates, which are useless.
    hData['X0'] = f.read('int16')
    hData['Y0'] = f.read('int16')
    # Number of data points alng x and y
    hData['Nx'] = f.read('int16')
    hData['Ny'] = f.read('int16')
    # Total data, 4*Nx*Ny
    hData['phaNBytes'] = f.read('int32')
    f.seek(218)
    # Scale factor determined by phase resolution tag
    phaseResTag = f.read('int16')
    if phaseResTag == 0:
        phaseResVal = 4096
    elif phaseResTag == 1:
        phaseResVal = 32768
    elif phaseResTag == 2:
        phaseResVal = 131072
    else:
        phaseResVal = 0
    f.seek(164)
    intf_scale_factor = f.read('float')
    hData['waveLength'] = f.read('float')
    f.seek(176)
    obliquity_factor = f.read('float')
    # Eq. in p12-6 in MetroPro Reference Guide
    hData['convFactor'] = intf_scale_factor*obliquity_factor*\
                          hData['waveLength']/phaseResVal
    # Bin size of each measurement
    f.seek(184)
    hData['cameraRes'] = f.read('float')
    
    return hData


class BinaryReaderEOFException(Exception):
    def __init__(self):
        pass
    def __str__(self):
        return 'Not enough bytes in file to satisfy read request'

class BinaryReader:
    # Map well-known type names into struct format characters.
    typeNames = {
        'int8'   :'b',
        'uint8'  :'B',
        'int16'  :'h',
        'uint16' :'H',
        'int32'  :'i',
        'uint32' :'I',
        'int64'  :'q',
        'uint64' :'Q',
        'float'  :'f',
        'double' :'d',
        'char'   :'s'}

    def __init__(self, fileName):
        self.file = open(fileName, 'rb')
        
    def read(self, typeName, size=None):
        import struct
        
        typeFormat = BinaryReader.typeNames[typeName.lower()]
        typeFormat = '>'+typeFormat
        typeSize = struct.calcsize(typeFormat)
        if size is None:
            value = self.file.read(typeSize)
            if typeSize != len(value):
                raise BinaryReaderEOFException
            unpacked = struct.unpack(typeFormat, value)[0]
        else:
            value = self.file.read(size*typeSize)
            if size*typeSize != len(value):
                raise BinaryReaderEOFException
            unpacked = np.zeros(size)
            for k in range(size):
                i = k*typeSize
                unpacked[k] = struct.unpack(typeFormat,value[i:i+typeSize])[0]
        return unpacked

    def seek(self, offset, refPos=0):
        '''
        offset in bytes and refPos gives reference position, where 0
        means origin of the file, 1 uses current position and 2 uses
        the end of the file.
        '''
        self.file.seek(offset, refPos)
    
    
    def __del__(self):
        self.file.close()


'''
TODO:

- map3=FT_invert_mirror_map(map3, invert) (Needed?)
    
- The intensity data is read but not used. Remove or add?

- When removing astigmatism and defocus separately, the process that is
  performed last will overwrite smap.Rc. Change this so that it is
  recognised that the other processes has been performed and calculate
  the "correct" Rc and store it.

- Similar to the above point, make zernikeRemoved additative so that
  when processing by fitting surfaces, and translating to Zernike,
  different processes can add to the same Zernike-polynomial amplitude.

- Messy with different index order for Zernike-polynomials. Tried to
  follow ddb's (m,n), but apparantly I have written (n,m) at a lot
  of places anyway. Nothing is wrong in the code though, buy a bit
  confusing when reading instructions and comments.

'''
    
