import os
import sys  
from flask import Flask
from optparse import OptionParser

repo_url = "https://git.ligo.org/finesse"

def start(instance_path,port=5000, debug=True, ip="0.0.0.0", git_bin="/usr/bin/git"):
        
    os.environ["GIT_BIN"] = git_bin
    # we import this now so that we can set the GIT_BIN env var
    from pykat.testing import utils
    
    print "starting web server..."
    
    if instance_path is None:
        raise Exception("instance_path must be defined")
    elif type(instance_path) is not str:
        raise Exception("instance_path must be a string")
        
    if not os.path.exists(instance_path):
        os.mkdir(instance_path)
        
    os.chdir(instance_path)
    
    from pykat.testing.web import app    
    
    if(app.instance_path!=instance_path):
	print app.instance_path, instance_path
        raise Exception("Instance path of Flask app (%s) didn't match the requested value (%s)" %(app.instance_path, instance_path))
    
    os.chdir(instance_path)    
    
    # need local copy of src
    if not os.path.exists(os.path.join(app.instance_path,"finesse_src")):
        print "finesse src folder didn't exist, cloning now..."
        utils.git(["clone","%s/finesse.git"%repo_url,"finesse_src"])
    else:
        # get the latest version for logs etc.
        utils.git("pull", cwd=os.path.join(app.instance_path,"finesse_src"))
     
    os.chdir(instance_path)
    
    # need local copy of test
    if not os.path.exists(os.path.join(app.instance_path,"finesse_test")):
        print "finesse test folder didn't exist, cloning now..."
        utils.git(["clone","%s/test.git"%repo_url,"finesse_test"])
        utils.git(["config","core.sharedRepository","true"], cwd="./finesse_test/")
    
    # load up the actual interface code
    import pykat.testing.web.web_interface
    
    app.secret_key = os.urandom(24)
    app.run(debug=debug, port=int(port), host=ip,use_reloader=False)

if __name__ == "__main__":
    
    parser = OptionParser()
    
    parser.add_option("-P","--path",type="string",dest="instance_path",help="")
    parser.add_option("-p","--port",type="int",default=5000,dest="port",help="")
    parser.add_option("-g","--git-bin",type="string",default="/usr/bin/git",dest="git_bin",help="")
    
    options, args = parser.parse_args()
    
    if options.instance_path is None:
        print "Must specify a path for the web server"
        exit()
    
    start(options.instance_path, port=options.port, git_bin=options.git_bin )
