# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 09:43:16 2013

@author: Daniel
"""

from distutils.core import setup

setup(
    name='PyKat',
    version='0.0.2',
    author='Daniel Brown',
    author_email='ddb@star.sr.bham.ac.uk',
    packages=['pykat','pykat.gui','pykat.gui.resources'],
    url='http://pypi.python.org/pypi/PyKat/',
    license='LICENSE.txt',
    description='Python interface and tools for FINESSE',
    long_description=open('README.txt').read(),
    install_requires=[
    "PyQt4 >= 4.8.3",
    "numpy >= 1.6.2",
    "flask >= 0.10.1"
    ],
)