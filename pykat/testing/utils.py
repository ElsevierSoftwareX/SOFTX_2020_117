import subprocess as sub
import os

GIT_BIN = "c:\\cygwin64\\bin\\git"

class RunException(Exception):
	def __init__(self, returncode, args, err, out):
		self.returncode = returncode
		self.args = args
		self.err = err
		self.out = out
        
def git(args):
    cmd = ""
    
    if type(args) is list:
        args.insert(0,GIT_BIN)
        cmd = " ".join(args)
    else:
        cmd = GIT_BIN + " " + args
        
    p = sub.Popen(cmd, stdout=sub.PIPE, stderr=sub.PIPE)
    out, err = p.communicate()

    if p.returncode != 0:
        print err
        raise RunException(p.returncode, args, err, out)

    return [out, err]
    
