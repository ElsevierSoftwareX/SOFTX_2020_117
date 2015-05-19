from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from pykat import finesse
from pykat.commands import *

def main():
		print("""
        --------------------------------------------------------------
        Template file for using PyKat to automate Finesse simulations
        Finesse: http://www.gwoptics.org/finesse
        PyKat:   http://www.gwoptics.org/pykat
                
        Andreas Freise 19.05.2015
        --------------------------------------------------------------
        """)   
        
        # for debugging we might need to see the temporay file:
		global kat
		kat = finesse.kat(tempdir=".",tempname="test")
		kat.verbose = False
		code1 = """
		l l1 1 0 0 n1
		s s1 1 n1 n2
		pd power n2
		xaxis l1 P lin 1 2 4
		"""
		kat.parseKatCode(code1)
		out = kat.run()
		print("result = {}".format(out.y))
	
if __name__ == '__main__':
        main()
