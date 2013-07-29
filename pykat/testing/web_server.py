import os
import sys  
from flask import Flask
from pykat.testing import utils

def start(instance_path,port=5000, debug=False):
    global app
    
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
        raise Exception("Instance path of Flask app didn't match the requested value")
        
    # load up the actual interface code
    import pykat.testing.web.web_interface
    
    # need local copy of src
    if not os.path.exists("./finesse_src"):
        print "finesse src folder didn't exist, cloning now..."
        utils.git("clone git://gitmaster.atlas.aei.uni-hannover.de/finesse/src.git finesse_src")
    else:
        os.chdir(os.path.join(app.instance_path,"finesse_src"))
        # get the latest version for logs etc.
        utils.git("pull")
        
    # need local copy of test
    if not os.path.exists("./finesse_test"):
        print "finesse test folder didn't exist, cloning now..."
        utils.git("clone git://gitmaster.atlas.aei.uni-hannover.de/finesse/test.git finesse_test")
    
    app.secret_key = os.urandom(24)
    app.run(debug=debug, port=int(port))

if __name__ == "__main__":
    n = len(sys.argv)
    
    if n != 2 and n != 3:
        print """
Command starts a Flask based website on which allows
users to manage the FINESSE test builds.
        
start_test_server usage:

    python -m pykat.testing.web_server <server_path> <port>

    Arguments:
        server_path: Directory where website should be run from
        port:        TCP port which the server should run under
                     default is port 5000.
"""
        exit()
        
    instance_path = sys.argv[1]
    
    if n == 3:
        port = sys.argv[2]
    else:
        port = 5000
    
    start(instance_path, port=port)
    