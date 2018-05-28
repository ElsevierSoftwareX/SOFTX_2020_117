#!/bin/bash

(cd .. && rm -r ./**/**.pyc)

(cd .. && python setup.py sdist --formats=zip upload)
