# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 09:43:16 2013

@author: Daniel
"""

from distutils.core import setup

REQUIREMENTS = [i.strip() for i in open("requirements.txt").readlines()]

setup(
    name='PyKat',
    version='0.2.4',
    author='Daniel Brown',
    author_email='ddb@star.sr.bham.ac.uk',
    packages=['pykat','pykat.gui','pykat.gui.resources','pykat.testing','pykat.testing.web'],
    url='http://pypi.python.org/pypi/PyKat/',
    license='GPL v2',
    description='Python interface and tools for FINESSE',
    long_description=open('README.rst').read(),
    install_requires=REQUIREMENTS
)