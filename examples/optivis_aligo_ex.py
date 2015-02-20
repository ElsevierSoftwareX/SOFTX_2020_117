import pykat

#kat = pykat.finesse.kat(kat_file="/Users/ddb/git/geo-simulation/geo-files/geoHF_v48.kat")

kat = pykat.finesse.kat(kat_file="LLO_matched.kat")

gui = kat.optivis()

gui.show()

