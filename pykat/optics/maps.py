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
from pykat.maths.zernike import *        
from pykat.exceptions import BasePyKatException

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
    def __init__(self, name, maptype, size, center, step_size, scaling, data=None,
                 notNan=None, Rc=None, zOffset=None):
        
        self.name = name
        self.type = maptype
        self.center = center
        self.step_size = step_size
        self.scaling = scaling
        self.notNan = notNan
        self.Rc = Rc
        self.zOffset = zOffset
        self.__interp = None
        
        if data is None:
            self.data = np.zeros(size)
        else:
            self.data = data
            
        if notNan is None:
            self.notNan = np.ones(size)
        else:
            self.notNan = notNan

        self._rom_weights = None
        
    def write_map(self, filename):
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
    def data(self):
        return self.__data
    
    @data.setter
    def data(self, value):
        self.__data = value
        self.__interp = None
    
    @property
    def center(self):
        return self.__center
    
    @center.setter
    def center(self, value):
        self.__center = value
        self.__interp = None
    
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
        return self.step_size[0] * (np.array(range(0, self.data.shape[0])) - self.center[0])
        
    @property
    def y(self):
        return self.step_size[1] * (np.array(range(0, self.data.shape[1]))- self.center[1])

    @property
    def size(self):
        return self.data.shape
            
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
    def plot(self, show=True, clabel=None, xlim=None, ylim=None):
        import pylab
        
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
            
        # ALSO ADDED (SEE LONG TEXT BELOW) BY DT TO FIX LIMITS
        # ------------------------------------------------------
        xlim,ylim = ylim,xlim
        # ------------------------------------------------------
        
        # min and max of z-values
        zmin = self.data[xmin:xmax,ymin:ymax].min()
        zmax = self.data[xmin:xmax,ymin:ymax].max()

        # 100 factor for scaling to cm
        xRange = 100*self.x
        yRange = 100*self.y
        
        # This line is added by DT to be able to plot
        # rectangular matrices. Effectively, I swapped the
        # x/y-axes. Preferrably, this should be corrected above
        # instead, but since I'm not completely sure of how the
        # coordinate system of these maps look I'll wait with
        # that. Here, I assume that row 0 of the matrix should
        # be plotted with y = Y[0], and that column 0 should be
        # plotted with x = X[0]. To be fully correct, I should
        # add one column and one row so that each matrix value
        # is plotted within the correct rectangle. 
        # ------------------------------------------------------
        xRange, yRange = np.meshgrid(yRange,xRange)
        # ------------------------------------------------------
        
        fig = pylab.figure()
        
        pcm = pylab.pcolormesh(xRange, yRange, self.data)
        pcm.set_rasterized(True)
        
        pylab.xlabel('x [cm]')
        pylab.ylabel('y [cm]')

        if xlim is not None: pylab.xlim(xlim)
        if ylim is not None: pylab.ylim(ylim)
            
        pylab.title('Surface map {0}, type {1}'.format(self.name, self.type))

        cbar = pylab.colorbar()
        #cbar.set_clim(zmin, zmax)
        
        if clabel is not None:
            cbar.set_label(clabel)
    
        if show:
            pylab.show()
        
        return fig



    def remove_curvature(self, Rc0, w=None, zOffset=None, isCenter=[False,False],
                         display='off'):
        '''
        Removes curvature from mirror map by fitting a sphere to
        mirror surface. Based on the file 'FT_remove_curvature_from_mirror_map.m'.
        Rc0      - Initial guess of the radius of curvature [m]
        w        - Beam radius on mirror [m], used for weighting.
        zOffset  - Initial guess of the z-offset [surfacemap.scaling]. Generally not needed.
        isCenter - 2D-list with booleans. isCenter[0] Determines if the center of the
                   sphere is fitted (True) or not (False, recommended). If the center is
                   fitted, then isCenter[1] determines if the weights are centered around
                   the fitted center (True) or centered around the center of the data-grid
                   (False, highly recommended).
        display  - Display mode of the fitting routine. Can be 'off',
                   'iter', 'notify', or 'final'.
        '''
        # z-value at centre of the data-grid. Serves as initial guess for deviation
        # from z=0, if no other first guess is given.
        if zOffset is None:
            zOffset = self.data[round(self.center[1]), round(self.center[0])]
            
        # If fitting center of the sphere, four variables are fitted. Initial guess
        # of deviation from notNan-data-grid-center: (x0,y0) = (0,0).
        if isCenter[0]:
            params = [Rc0, zOffset, 0.0, 0.0]
        # Otherwise two.
        else:
            params = [Rc0, zOffset]

        # Grid with X,Y and r2 = (X^2+Y^2) values. X and Y crosses zero in the center
        # of the xy-plane.
        X,Y,r2 = self.createGrid()
        
        # Cost-function to minimize.
        def costFunc(params):
            # If the center of the mirror is fitted, four variables are fitted.
            if isCenter[0]:
                Rc = params[0]
                zOffset = params[1]
                x0 = params[2]
                y0 = params[3]
            # Otherwise 2 variables.
            else:
                Rc = params[0]
                zOffset = params[1]
                x0 = 0
                y0 = 0
            
            Z = self.createSphere(Rc,X,Y,zOffset,x0,y0)
            
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
                res = math.sqrt( (weight[self.notNan]*((self.data[self.notNan] -
                                   Z[self.notNan])**2)).sum() )/weight[self.notNan].sum()
            return res
        
        # Using the simplex Nelder-Mead method. This is the same or very
        # similar to the method used in 'FT_remove_curvature_from_mirror_map.m',
        # but there are probably better methods to use.
        out = minimize(costFunc, params, method='Nelder-Mead',tol=1.0e-6)

        # Assigning values to the instance variables
        self.Rc = out['x'][0]
        self.zOffset = out['x'][1]

        # If center was fitted, assign new values to instance variable center, and
        # subtract the fitted sphere from the mirror map.
        if isCenter[0]:
            x0 = out['x'][2]
            y0 = out['x'][3]
            # Converts the deviation into a new absolut center in data points.
            self.center = (self.center[0] + x0/self.step_size[0],
                           self.center[1] + y0/self.step_size[1])
            # Creating fitted sphere
            Z = self.createSphere(self.Rc, X, Y, self.zOffset, x0, y0)
            # Subtracting sphere from map
            self.data[self.notNan] = self.data[self.notNan]-Z[self.notNan]
            return self.Rc, self.zOffset, x0, y0
        # Subtracting fitted sphere from mirror map.
        else:
            # Creating fitted sphere
            Z = self.createSphere(self.Rc,X,Y,self.zOffset)
            # Subtracting sphere from map
            self.data[self.notNan] = self.data[self.notNan]-Z[self.notNan]
            return self.Rc, self.zOffset

    def createSphere(self,Rc,X,Y,zOffset=0,x0=0,y0=0,xTilt=0,yTilt=0,isPlot=False):
        '''
        Creating spherical surface. Based on 'FT_create_sphere_for_map.m'
        Rc      - Radius of curvature, and center of sphere on z-axis in case zOffset=0 [m].
        X, Y    - Matrices of x,y-values generated by 'X,Y = numpy.meshgrid(xVec,yVec)'.
                  Preferrably created by method 'surfacemap.createGrid()' [m].
        zOffset - Surface center offset from 0. Positive value gives the surface center
                  positive z-coordinate. [self.scaling]
        x0,y0   - The center of the sphere in the xy-plane.
        x/yTilt - Tilt of x/y-axis [rad].
        isPlot  - Plot or not [boolean]. Not recommended when used within optimization
                  algorithm.
        '''
        
        # Adjusting for tilts and offset
        Z = zOffset + (X*np.tan(xTilt) + Y*np.tan(yTilt))/self.scaling
        # Adjusting for spherical shape.
        Z = Z + (Rc - np.sqrt(Rc**2 - (X-x0)**2 - (Y-y0)**2))/self.scaling

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
    
    def createGrid(self):
        '''
        Creating grid for fitting spherical surface to the map. Based on
        'FT_create_grid_for_map.m' and 'FT_init_grid.m'. (x,y) = (0,0) is
        placed at the centre of the mirror surface defined by the instance
        variable surfacemap.center.
        
        Returns X,Y,r2
        , where X,Y = numpy.meshgrid(x,y), and r2 = X^2 + Y^2.
        '''
        
        yPoints, xPoints = self.data.shape
        # Difference between the data point centre, and the centre of
        # the mirror surface.
        xOffset = (((xPoints-1)/2-self.center[0]) )*self.step_size[0]
        yOffset = (((yPoints-1)/2-self.center[1]) )*self.step_size[1]

        # Creating arrays with values from 0 up to (xyPoints-1).
        x = np.array([n for n in range(xPoints)])
        y = np.array([n for n in range(yPoints)])

        # Physical size. Shouldn't this be (points-1)*step_size?
        xSize = xPoints*self.step_size[0]
        ySize = yPoints*self.step_size[1]
        # ...corrected here by subracting one step. Isn't this unnecessary
        # complicated btw?
        xAxis = x*self.step_size[0] + xOffset - (xSize-self.step_size[0])/2
        yAxis = y*self.step_size[1] + yOffset - (ySize-self.step_size[1])/2
        
        X,Y = np.meshgrid(xAxis,yAxis)
        r2 = X**2 + Y**2
        # r = np.sqrt(X**2 + Y**2)
        # phi = np.arctan2(Y,X)
        return X,Y,r2
        # ----------------------------------------------------------------

        

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

    def plot(self, mode="absorption", show=True, clabel=None, xlim=None, ylim=None, wavelength=1064e-9):
        
        import pylab
        
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

        fig = pylab.figure()
        axes = pylab.pcolormesh(xrange, yrange, data, vmin=zmin, vmax=zmax)
        pylab.xlabel('x [cm]')
        pylab.ylabel('y [cm]')

        if xlim is not None: pylab.xlim(xlim)
        if ylim is not None: pylab.ylim(ylim)

        pylab.title('Merged map {0}, mode {1}'.format(self.name, mode))

        cbar = fig.colorbar(axes)
        cbar.set_clim(zmin, zmax)
        
        if clabel is not None:
            cbar.set_label(clabel)
    
        if show:
            pylab.show()
        
        return fig

class aperturemap(surfacemap):
    
    def __init__(self, name, size, step_size, R):
        surfacemap.__init__(self, name, "absorption both", size, (np.array(size)+1)/2.0, step_size, 1)
        self.R = R
        
    @property
    def R(self):
        return self.__R
    
    @R.setter
    def R(self, value):
        self.__R = value
    
        xx, yy = np.meshgrid(self.x, self.y)
        
        radius = np.sqrt(xx**2 + yy**2)
        
        self.data = np.zeros(self.size)
        self.data[radius > self.R] = 1.0
        
        
class curvedmap(surfacemap):
    
    def __init__(self, name, size, step_size, Rc):
        surfacemap.__init__(self, name, "phase reflection", size, (np.array(size)+1)/2.0, step_size, 1e-6)
        self.Rc = Rc
        
    @property
    def Rc(self):
        return self.__Rc
    
    @Rc.setter
    def Rc(self, value):
        self.__Rc = value
    
        xx, yy = np.meshgrid(self.x, self.y)
        
        Rsq = xx**2 + yy**2
        self.data = (self.Rc - math.copysign(1.0, self.Rc) * np.sqrt(self.Rc**2 - Rsq))/ self.scaling

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
	
			

def read_map(filename, mapFormat='finesse', scaling=1.0e-9):
    '''
    Reads surface map files and return a surfacemap-object xy-centered at the
    center of the mirror surface.
    filename  - name of surface map file.
    mapFormat - 'finesse', 'ligo', 'zygo'. Currently only for ascii formats.
    scaling   - scaling of surface height of the mirror map [m].
    '''
    
    # Function converting string number into float.
    g = lambda x: float(x)

    # Reads finesse mirror maps.
    if mapFormat == 'finesse':
        
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
        
    # Converts raw zygo and ligo mirror maps to the finesse
    # format. Based on the matlab scripts 'FT_read_zygo_map.m' and
    # 'FT_read_ligo_map.m'.
    elif mapFormat == 'ligo' or mapFormat == 'zygo':
        if mapFormat == 'ligo':
            isLigo = True
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
                
        # Unknowns (why are these values hard coded here?) Moving
        # scaling to input. Fix the others later too.
        # ------------------------------------------------------
        # Standard maps have type 'phase' (they store surface
        # heights)
        maptype = 0
        # Both (reflected and transmitted) light fields are
        # affected
        field = 0
        # Measurements in nanometers
        # scaling = 1.0e-9
        # ------------------------------------------------------

        # Reading header of LIGO-map (Zygo file? Says Zygo in
        # header...)
        # ------------------------------------------------------
        with open(filename, 'r') as f:
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
            # Unknown usage.
            # ----------------------------------------------
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
            # ----------------------------------------------

            # Skipping three lines
            for k in range(3):
                f.readline()
            line = f.readline().split()

            # Unknown (Scaling factors)
            # ----------------------------------------------
            # Interfeometric scaling factor (?)
            S = float(line[1])
            # wavelength (of what?)
            lam = float(line[2])
            # Obliquity factor (?)
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

            # Unknown
            # Resolution of phase data points, 1 or 0.
            phaseRes = float(line[0])
            if phaseRes == 0:
                R = 4096
            elif phaseRes == 1:
                R = 32768
            else:
                print('Error, invalid phaseRes')

            if not isLigo and not isAscii:
                # zygo .xyz files give phase data in microns.
                hScale = 1.0e-6
            else:
                # zygo .asc and ligo-files give phase data in
                # internal units. To convert to m use hScale
                # factor.
                hScale = S*O*lam/R
                
            if not isLigo and not isAscii:
                print('Not implemented yet, need a .xyz-file ' +
                      'to do this.')
                return 0
                
            # Skipping four lines
            for k in range(4):
                f.readline()
            if not isLigo and isAscii:
                # Reading intensity data
                iData = np.array([])
                line = f.readline().split()
                while line[0] != '#':
                    iData = np.append(iData, map(g,line))
                    line = f.readline().split()
                # Reshaping intensity data
                iData = iData.reshape(iRows, iCols).transpose()
                iData = np.rot90(iData)
            else:
                # Skipping lines until '#' is found.
                while f.readline()[0] != '#':
                    pass
                
            # Reading phase data
            # ----------------------------------------------
            # Array with the data
            data = np.array([])
            # Reading data until next '#' is reached.
            line = f.readline().split()
            while line[0] != '#':
                data = np.append(data, map(g,line))
                line = f.readline().split()
            # ----------------------------------------------

        
        if isLigo:
            # Setting all the points outside of the mirror
            # surface to NaN. These are given a large number
            # in the file. 
            data[data == data[0]] = np.nan
            
            # Reshaping into rows and columns
            data = data.reshape(cols,rows).transpose()
            # Pretty sure that the lines below can be done
            # more efficient, but it's quick as it is.
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
            data = data.reshape(rows,cols).transpose()
            # Rotating to make (0,0) be in bottom left
            # corner. 
            data = np.rot90(data)
            
        # Scaling to nanometer (change this to a user
        # defined value?) Still don't know where
        # 'hScale' really comes from.
        data = (hScale/scaling)*data
        size = data.shape

        if maptype == 0:
            mType = 'phase'
        else:
            mType = 'Unknown'
        if field == 0:
            fType = 'both'
        else:
            fType = 'unknown'

        maptype = ' '.join([mType, fType])

        # Wrong! fix by creating recenter method.
        center = tuple([x0,y0])
        step = tuple([xstep,ystep])

        # Simple re-centering of mirror, translated from
        # 'FT_recenter_mirror_map.m'
        # -------------------------------------------------
        # Matrix with ones where data element is not NaN.
        isNan = np.isnan(data)
        notNan = isNan==False
        # Row and column indices with non-NaN elements
        rIndx, cIndx = notNan.nonzero()
        # Finding centres
        x0 = float(cIndx.sum())/len(cIndx)
        y0 = float(rIndx.sum())/len(rIndx)
        center = tuple([x0,y0])
        # -------------------------------------------------
        
        # Changing NaN to zeros. Just to be able to plot the
        # map with surfacemap.plot().
        data[isNan] = 0 
        
        return surfacemap(name, maptype, size, center, step, scaling, data, notNan)
        
    # TODO: Add options for reading virgo maps, and .xyz zygo
    # maps (need .xys file for this). Binary ligo-maps?
    # The intensity data is not used to anything here. Remove
    # or add to pykat?

    
    


# TODO: Recreate functions from Simtools:, List taken from: ligo_maps/FT_convert_ligo_map_for_finesse.m
# map=FT_recenter_mirror_map(map);
# [map2,A2,Rc_out]=FT_remove_zernike_curvatures_from_map(map,Rc_in);
# [map2,Rc_out]=FT_remove_curvature_from_mirror_map(map,Rc_in,w, display_style);
# [map2,offset]=FT_remove_offset_from_mirror_map(map2,1e-2);
# [map3,x_tilt,y_tilt,offset2]=FT_remove_piston_from_mirror_map(map2,w, display_style);
# map3=FT_invert_mirror_map(map3, invert);

# Understand the internal coordinate system of the
# maps/matrices.
