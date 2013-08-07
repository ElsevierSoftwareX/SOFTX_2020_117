import subprocess as sub
import os

try:
    GIT_BIN = os.environ["GIT_BIN"]
except:
    GIT_BIN = "/usr/bin/git"

print "GIT_BIN = " + GIT_BIN
    
class RunException(Exception):
	def __init__(self, returncode, args, err, out):
		self.returncode = returncode
		self.args = args
		self.err = err
		self.out = out

def runcmd(args, cwd="."):
    p = sub.Popen(args, stdout=sub.PIPE, stderr=sub.PIPE, cwd=cwd)
    out, err = p.communicate()
    
    if p.returncode != 0:
        print "STDERR: " + err
        print "STDOUT: " + out
        raise RunException(p.returncode, args, err, out)

    return [out,err]
    
def git(args, git_bin=GIT_BIN, cwd="."):
    cmd = ""

    if type(args) is not list:
        if type(args) is not str:
            print type(args)
            raise Exception("arg for utils.git must be a list or string")
    
    if type(args) is str:
        args = args.split(" ")
    
    args.insert(0, git_bin)
    
    p = sub.Popen(args, stdout=sub.PIPE, stderr=sub.PIPE, cwd=cwd)
    out, err = p.communicate()
        
    if p.poll() != 0:
        raise RunException(p.poll(), args, err, out)

    return [out, err]
    
