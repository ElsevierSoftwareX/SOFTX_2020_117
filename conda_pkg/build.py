import os
import os.path
import subprocess



subprocess.check_call(["conda", "config", "--add", "channels", "gwoptics"])
subprocess.check_call(["conda", "build", "purge"])
subprocess.check_call(["conda", "build", "--output-folder", "packages", "pykat"])
    
    