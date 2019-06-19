#!/bin/bash

(cd .. && rm -r ./**/**.pyc)

(cd .. && python setup.py sdist --formats=zip)

(cd .. && twine upload dist/*)

