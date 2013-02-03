# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 09:43:16 2013

@author: Daniel
"""

from distutils.core import setup

setup(
    name='PyKat',
    version='0.0.1',
    author='Daniel Brown',
    author_email='danielbrown87@gmail.com',
    packages=['pykat'],
    url='http://www.gwoptics.org/pykat',
    license='LICENSE.txt',
    description='Python interface and tools for FINESSE',
    long_description=open('README.txt').read(),
    install_requires=[
    ],
)