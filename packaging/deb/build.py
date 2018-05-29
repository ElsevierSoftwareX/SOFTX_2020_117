import os
from subprocess import call, check_output

class SafeDict(dict):
    def __missing__(self, key):
        return '{' + key + '}'

os.chdir("/root/pykat")

call("git pull".split())

git_describe = str(check_output('git describe --tags'.split())).split("-")

vals = {
    "version": git_describe[0],
    "release": git_describe[1],
}

print(vals)

os.chdir("/root")            
call("fpm -s python -t deb pykat/setup.py".split())
call("cp python-pykat_{version}.{release}_all.deb /host".format(**vals).split())