from pykat.utilities.romhom import makeReducedBasis, makeEmpiricalInterpolant, makeWeights
import numpy
import math

class surfacemap:
    def __init__(self, name, maptype, size, center, step_size, scaling, data=None):
        
        self.name = name
        self.type = maptype
        self.center = center
        self.step_size = step_size
        self.scaling = scaling
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
    def x(self):
        return self.step_size[0] * (numpy.array(range(1, self.data.shape[0]+1)) - self.center[0])
        
    @property
    def y(self):
        return self.step_size[1] * (numpy.array(range(1, self.data.shape[1]+1))- self.center[1])

    @property
    def size(self):
        return self.data.shape
            
    @property
    def offset(self):
        return numpy.array(self.step_size)*(self.center - numpy.array(self.size)/2)
    
    @property
    def ROMWeights(self):
        return self._rom_weights
    
    def z_xy(self, wavelength=1064e-9):
        
        if "phase" in self.type:
            k = math.pi * 2 / wavelength
            return numpy.exp(2j * k * self.scaling * self.data)
        else:
            raise BasePyKatException("Map type needs handling")
        
    
    def generateROMWeights(self):
        b = makeReducedBasis(self.x[0:(len(self.x)/2)], offset=self.offset)
        EI = makeEmpiricalInterpolant(b)
        self._rom_weights = makeWeights(self, EI)
        
        return self.ROMWeights, EI
        
    def plot(self, show=True, clabel=None):
        
        import pylab
        
        # 100 factor for scaling to cm
        xrange = 100*self.x
        yrange = 100*self.y
        
        fig = pylab.figure()
        axes = pylab.imshow(self.data, extent=[min(xrange),max(xrange),min(yrange),max(yrange)])
        pylab.xlabel('x [cm]')
        pylab.ylabel('y [cm]')

        pylab.title('Surface map {0}, type {1}'.format(self.name, self.type))
        
        cbar = fig.colorbar(axes)
        
        if clabel != None:
            cbar.set_label(clabel)
                
        if show:
            pylab.show()
            
        return fig

        
  
def read_map(filename):
    with open(filename, 'r') as f:
        
        f.readline()
        name = f.readline().split(':')[1].strip()
        maptype = f.readline().split(':')[1].strip()
        size = tuple(map(lambda x: int(x), f.readline().split(':')[1].strip().split()))
        center = tuple(map(lambda x: float(x), f.readline().split(':')[1].strip().split()))
        step = tuple(map(lambda x: float(x), f.readline().split(':')[1].strip().split()))
        scaling = float(f.readline().split(':')[1].strip())
        
        
        
    data = numpy.loadtxt(filename, dtype=numpy.float64,ndmin=2,comments='%')    
        
    return surfacemap(name,maptype,size,center,step,scaling,data)
    
    
        
        