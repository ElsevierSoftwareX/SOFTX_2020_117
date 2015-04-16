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

from pykat.optics.romhom import makeReducedBasis, makeEmpiricalInterpolant, makeWeights, makeWeightsNew
from scipy.interpolate import interp2d, interp1d
from pykat.maths.zernike import *        

import numpy as np
import math
import pickle
		   
class surfacemap(object):
    def __init__(self, name, maptype, size, center, step_size, scaling, data=None):
        
        self.name = name
        self.type = maptype
        self.center = center
        self.step_size = step_size
        self.scaling = scaling
        self.__interp = None
        
        if data is None:
            self.data = np.zeros(size)
        else:
            self.data = data

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
        return self.step_size[0] * (np.array(range(1, self.data.shape[0]+1)) - self.center[0])
        
    @property
    def y(self):
        return self.step_size[1] * (np.array(range(1, self.data.shape[1]+1))- self.center[1])

    @property
    def size(self):
        return self.data.shape
            
    @property
    def offset(self):
        return np.array(self.step_size)*(np.array(self.center) - 1/2. - np.array(self.size)/2.0)
    
    @property
    def ROMWeights(self):
        return self._rom_weights
    
    def z_xy(self, x=None, y=None, wavelength=1064e-9, direction="reflection", nr1=1.0, nr2=1.0):
        
        if x is None and y is None:
            data = self.scaling * self.data
        else:
            if self.__interp is None:
                self.__interp = interp2d(self.x, self.y, self.data * self.scaling)
                
            data = self.__interp(x, y)
            
        if direction == "reflection":
            if "phase" in self.type:
                k = math.pi * 2 / wavelength
                return np.exp(2j * k * data)
                
            elif "absorption" in self.type:
                return np.sqrt(1.0 - data)
            else:
                raise BasePyKatException("Map type needs handling")
                
        elif direction == "transmission":
            if "phase" in self.type:
                k = math.pi * 2 / wavelength
                return np.exp((nr1-nr2)*k * data)
                
            elif "absorption" in self.type:
                return np.sqrt(1.0 - data)
                
            else:
                raise BasePyKatException("Map type needs handling")
        else:
            raise BasePyKatException("Map type needs handling")
        

    def generateROMWeights(self, isModeMatched=True, verbose=False, interpolate=False, interpolate_N=None, tolerance = 1e-12, sigma = 1, sort=False, greedyfile=None, useSymmetry=True):
        
        if interpolate:
            from scipy.interpolate import interp2d
            import numpy as np

            D = interp2d(self.x, self.y, self.data, fill_value=0)
            if interpolate_N is None:
                interpolate_N = (self.size[0], self.size[1])
                
            # only want even number of data points spread equally
            # about the axes
            if interpolate_N[0] % 2 == 0:
                Nx = interpolate_N[0]
            else:
                Nx = interpolate_N[0]-1

            if interpolate_N[1] % 2 == 0:
                Ny = interpolate_N[1]
            else:
                Ny = interpolate_N[1]-1
            
            nx = np.linspace(min(self.x), max(self.x), interpolate_N[0]) 
            ny = np.linspace(min(self.y), max(self.y), interpolate_N[1])
            
            data = D(nx-self.offset[0], ny-self.offset[0])
            
            self.name += "[ROMHOM_Interpolated]"
            
            self.center = (np.array(data.shape)+1)/2.0
            
            self.data = data
    
        if useSymmetry:
            xm = self.x[self.x<0]
            ym = self.y[self.y<0]
        else:
            xm = self.x
            ym = self.y
        
        if min(xm) == min(ym) and max(xm) == max(ym) and len(xm) == len(ym):
            symm = True
        else:
            symm = False
            
        EI = {}
        
        EI["x"] = makeEmpiricalInterpolant(makeReducedBasis(xm, isModeMatched=isModeMatched, tolerance = tolerance, sigma = sigma, greedyfile=greedyfile), sort=sort)
        
        if symm:
            EI["y"] = EI["x"]
        else:
            EI["y"] = makeEmpiricalInterpolant(makeReducedBasis(ym, isModeMatched=isModeMatched, tolerance = tolerance, sigma = sigma, greedyfile=greedyfile), sort=sort)
        
        EI["limits"] = EI["x"].limits
        
        self._rom_weights = makeWeights(self, EI, verbose=verbose, useSymmetry=useSymmetry)
        
        return self.ROMWeights, EI
    
    
    def generateROMWeightsNew(self, EIxFilename, EIyFilename=None, verbose=False, interpolate=False):
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
        
        self._rom_weights = makeWeightsNew(self, EIxFilename, EIyFilename, verbose=verbose)
        return self.ROMWeights

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
        
        data = D(nx-self.offset[0], ny-self.offset[0])
        
        Dx = interp1d(nx, np.arange(1,len(nx)+1))
        Dy = interp1d(ny, np.arange(1,len(ny)+1))
        
        self.center = (Dx(0), Dy(0))
        self.step_size = (nx[1]-nx[0], ny[1]-ny[0])
        self.data = data

    def plot(self, show=True, clabel=None, xlim=None, ylim=None):
        
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
        
        zmin = self.data[xmin:xmax,ymin:ymax].min()
        zmax = self.data[xmin:xmax,ymin:ymax].max()

        # 100 factor for scaling to cm
        xrange = 100*self.x
        yrange = 100*self.y

        fig = pylab.figure()
        axes = pylab.pcolormesh(xrange, yrange, self.data, vmin=zmin, vmax=zmax)
        pylab.xlabel('x [cm]')
        pylab.ylabel('y [cm]')

        if xlim is not None: pylab.xlim(xlim)
        if ylim is not None: pylab.ylim(ylim)

        pylab.title('Surface map {0}, type {1}'.format(self.name, self.type))

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
	
			
	
def read_map(filename):
    with open(filename, 'r') as f:
        
        f.readline()
        name = f.readline().split(':')[1].strip()
        maptype = f.readline().split(':')[1].strip()
        size = tuple(map(lambda x: int(x), f.readline().split(':')[1].strip().split()))
        center = tuple(map(lambda x: float(x), f.readline().split(':')[1].strip().split()))
        step = tuple(map(lambda x: float(x), f.readline().split(':')[1].strip().split()))
        scaling = float(f.readline().split(':')[1].strip())
        
        
        
    data = np.loadtxt(filename, dtype=np.float64,ndmin=2,comments='%')    
        
    return surfacemap(name,maptype,size,center,step,scaling,data)
    
    
        
        
