import os
import os.path
import subprocess
import pykat

import subprocess

version_git_1 = subprocess.check_output(["git", "describe","--long"]).decode('utf8').rstrip()
version_git = ".".join(version_git_1.split('-')[:2])
    
os.environ['CONDA_PYKAT_VERSION'] = version_git
os.environ['CONDA_BUILD_VERSION'] = version_git_1

subprocess.check_call(["conda", "config", "--add", "channels", "gwoptics"])
subprocess.check_call(["conda", "build", "purge"])
subprocess.check_call(["conda", "build", "--output-folder", "packages", "pykat"])
    
    