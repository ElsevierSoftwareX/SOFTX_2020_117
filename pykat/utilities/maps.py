from pykat.utilities.romhom import makeReducedBasis, makeEmpiricalInterpolant, makeWeights
from scipy.interpolate import interp2d
import numpy as np
import math
        
        
class surfacemap(object):
    def __init__(self, name, maptype, size, center, step_size, scaling, data=None):
        
        self.name = name
        self.type = maptype
        self.center = center
        self.step_size = step_size
        self.scaling = scaling
        self.__interp = None
        
        if data == None:
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
        
        if x == None and y == None:
            data = self.scaling * self.data
        else:
            if self.__interp == None:
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
        

    def generateROMWeights(self, isModeMatched=True, verbose=False, interpolate=False, tolerance = 1e-12, sigma = 1, sort=False, greedyfile=None):
        
        if interpolate:
            from scipy.interpolate import interp2d
            import numpy as np

            D = interp2d(self.x, self.y, self.data, fill_value=0)

            # only want even number of data points spread equally
            # about the axes
            if self.size[0] % 2 == 0:
                Nx = self.size[0]
            else:
                Nx = self.size[0]-1

            if self.size[1] % 2 == 0:
                Ny = self.size[1]
            else:
                Ny = self.size[1]-1
            
            nx = np.linspace(min(self.x), max(self.x), Nx) 
            ny = np.linspace(min(self.y), max(self.y), Ny)
            
            data = D(nx-self.offset[0], ny-self.offset[0])
            
            self.name += " [ROMHOM interpolated]"
            
            self.center = (np.array(data.shape)+1)/2.0
            
            self.data = data
        
        xm = self.x[self.x<0]
        ym = self.y[self.y<0]
        
        symm = False
        
        if min(xm) == min(ym) and max(xm) == max(ym) and len(xm) == len(ym):
            symm = True
        
        EI = {}
        
        EI["xm"] = makeEmpiricalInterpolant(makeReducedBasis(xm, isModeMatched=isModeMatched, tolerance = tolerance, sigma = sigma, greedyfile=greedyfile), sort=sort)
        
        if symm:
            EI["ym"] = EI["xm"]
        else:
            EI["ym"] = makeEmpiricalInterpolant(makeReducedBasis(ym, isModeMatched=isModeMatched, tolerance = tolerance, sigma = sigma, greedyfile=greedyfile), sort=sort)
        
        EI["limits"] = EI["xm"].limits
        
        self._rom_weights = makeWeights(self, EI, verbose=verbose)
        
        return self.ROMWeights, EI
        
    def plot(self, show=True, clabel=None, xlim=None, ylim=None):
        
        import pylab
        
        if xlim != None:
            _x = np.logical_and(self.x<=max(xlim)/100.0, self.x>=min(xlim)/100.0)
            xmin = np.min(np.where(_x == True))
            xmax = np.max(np.where(_x == True))
        else:
            xmin = 0
            xmax = len(self.x)-1
            xlim = [self.x.min()*100, self.x.max()*100]
    
        if ylim != None:
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

        if xlim != None: pylab.xlim(xlim)
        if ylim != None: pylab.ylim(ylim)

        pylab.title('Surface map {0}, type {1}'.format(self.name, self.type))

        cbar = fig.colorbar(axes)
        cbar.set_clim(zmin, zmax)
        
        if clabel != None:
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
        surfacemap.__init__(self, name, "phase", size, (np.array(size)+1)/2.0, step_size, 1e-9)
        self.tilt = tilt
        
    @property
    def tilt(self):
        return self.__tilt
    
    @tilt.setter
    def tilt(self, value):
        self.__tilt = value
        
        xx, yy = np.meshgrid(self.x, self.y)
        
        self.data = (xx * self.tilt[1] + yy * self.tilt[0])/self.scaling
        
        
        
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
    
    
        
        