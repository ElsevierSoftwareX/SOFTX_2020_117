# -*- coding: utf-8 -*-
# Author: Aaron Jones
# To test KatRun2D save and load
# Simple coupled cavity as test case

import pykat
from pykat import finesse

def test_string(_str):
    kat = finesse.kat() # create a fresh cat object
    _str = _str.strip()
    kat.parse(_str)

    # remove extra kat info produced by pykat
    katfile = ''.join(kat.generateKatScript()[1:-1]).strip()
    try:
        assert katfile == _str
    except AssertionError as e:
        print('Input string {} does not match output string {}'.format(_str,katfile))
        raise e

# Check that lens are read and output correctly
test_string("lens l1 10.0 n0 n1")
test_string("lens* l1 10.0 n0 n1")
test_string("lens** l1 10.0 15.0 n0 n1")
test_string("lens*** l1 10.0 15.0 n0 n1")

# Check the pykat interface works
kat = finesse.kat() # create a fresh cat object
kat.parse("lens l1 10.0 n0 n1")
assert kat.l1.p == None
assert kat.l1.f == 10.0

kat = finesse.kat() # create a fresh cat object
kat.parse("lens* l1 10.0 n0 n1")

assert kat.l1.p == 10.0
assert kat.l1.f == None

kat = finesse.kat() # create a fresh cat object
kat.parse("lens** l1 10.0 15.0 n0 n1")

assert kat.l1.px == None
assert kat.l1.py == None
assert kat.l1.fx == 10.0
assert kat.l1.fy == 15.0

kat = finesse.kat() # create a fresh cat object
kat.parse("lens*** l1 10.0 15.0 n0 n1")

assert kat.l1.fx == None
assert kat.l1.fy == None
assert kat.l1.px == 10.0
assert kat.l1.py == 15.0

