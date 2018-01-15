import pykat
from pykat.optics.maps import read_map
from pykat.optics.maps import curvedmap

m = read_map("etm08_virtual.txt")

m.write_map("test.map")

itm = curvedmap('itm_Rc', (10,10), (1,1), 100)



