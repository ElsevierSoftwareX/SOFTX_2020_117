import os
import os.path
import subprocess

os.environ['CONDA_BUILD_VERSION'] = 'master'

subprocess.check_call(["conda", "config", "--add", "channels", "gwoptics"])

for _ in ["2.7", "3.4", "3.5", "3.6"]:
    subprocess.check_call(["conda", "build", "purge"])
    subprocess.check_call(["conda", "build", "--py", _, "pykat"])