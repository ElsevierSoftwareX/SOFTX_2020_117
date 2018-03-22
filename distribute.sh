#!/bin/bash

rm -r ./**/**.pyc

python setup.py sdist --formats=gztar,zip upload
