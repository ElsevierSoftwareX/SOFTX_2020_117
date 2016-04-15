# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 09:43:16 2013

@author: Daniel
"""
from pykat import __version__ as version
from distutils.core import setup
import os

REQUIREMENTS = [i.strip() for i in open("requirements.txt").readlines()]

setup(
    name='PyKat',
    version=version,
    author='Daniel Brown',
    author_email='ddb@star.sr.bham.ac.uk',
    packages=[x[0].replace("/",".") for x in os.walk("pykat") if "__" not in x[0]],
    url='http://pypi.python.org/pypi/PyKat/',
    license='GPL v2',
    description='Python interface and tools for FINESSE',
    long_description=open('README.rst').read(),
    install_requires=REQUIREMENTS,
    package_data={'': ['optics/greedypoints/*.txt']},
    include_package_data=True
)
