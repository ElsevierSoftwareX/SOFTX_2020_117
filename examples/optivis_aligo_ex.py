import pykat

kat = pykat.finesse.kat(kat_file="/Users/ddb/git/geo-simulation/geo-files/geoHF_v48.kat")

gui, scene = kat.optivis()

gui.show()

