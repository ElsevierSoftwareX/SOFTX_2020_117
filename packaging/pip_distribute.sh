#!/bin/bash

(cd .. && rm -r ./**/**.pyc)

(cd .. && python setup.py sdist --formats=gztar,zip)

(cd .. && twine upload dist/*)

