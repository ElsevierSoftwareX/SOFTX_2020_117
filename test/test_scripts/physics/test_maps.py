import pykat
from pykat.optics.maps import read_map

m = read_map("etm08_virtual.txt")

m.write_map("test.map")



from pykat.optics.maps import curvedmap

itm = curvedmap('itm_Rc', 1, 1, 1)



