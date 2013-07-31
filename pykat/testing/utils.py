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

def runcmd(args):
    p = sub.Popen(args, stdout=sub.PIPE, stderr=sub.PIPE)
    out, err = p.communicate()
    
    if p.returncode != 0:
        print "STDERR: " + err
        print "STDOUT: " + out
        raise RunException(p.returncode, args, err, out)

    return [out,err]
    
def git(args, git_bin=GIT_BIN):
    cmd = ""
    
    if type(args) is list:
        args.insert(0,git_bin)
        cmd = " ".join(args)
    else:
        cmd = git_bin + " " + args
        
    print cmd
    
    print os.getcwd()
    
    p = sub.Popen(cmd, stdout=sub.PIPE, stderr=sub.PIPE)
    out, err = p.communicate()
        
    if p.returncode != 0:
        print "STDERR: " + err
        print "STDOUT: " + err
        raise RunException(p.poll(), args, err, out)

    return [out, err]
    
