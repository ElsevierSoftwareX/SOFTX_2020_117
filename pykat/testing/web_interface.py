from threading import Thread, Lock
from time import sleep
from uuid import uuid4
from flask import Flask
from flask import jsonify
from flask import Module
from flask import request
from flask import render_template
from datetime import datetime
from collections import namedtuple

import utils
import os

global current_test, scheduled_tests, schedule_lock

test_id = 0
current_test = None
scheduled_tests = []
schedule_lock = Lock()

app = Flask(__name__)

@app.route('/')
def hello_world():
    return render_template("finesse_test.html")

@app.route('/finesse/start_test', methods=["POST"])
def finesse_start_test():
    global current_test, test_id
    
    try:
        schedule_lock.acquire()
        
        test_id += 1
        
        git_commit = request.json['git_commit']
                
        test = FinesseTestProcess(git_commit, test_id)
        
        # check if anything is running and if it
        # isn't start this test off
        if current_test is None:
            print "running test"
            current_test = test
            
            test.start()
        else:
            print "queuing test"
            scheduled_tests.append(test)
    finally:
        schedule_lock.release()
    
    return jsonify({'id':test.test_id, 'queued': (current_test != test)})
        
@app.route('/finesse/get_tests', methods=["POST"])
def finesse_get_tests():
    tests = []
    
    for test in scheduled_tests:
        tests.append({'id':test.test_id,'git_commit':test.get_version(),'queue_time':test.queue_time});
    
    return jsonify(tests=tests)
    
@app.route('/finesse/get_test_progress', methods=["POST"])
def finesse_get_test_progress():
    percent_done = 0
    status = ""
    version = ""
    
    try:
        schedule_lock.acquire()
    
        if current_test is None:
    
            return jsonify(running=False)
        else:
    
            percent_done = current_test.percent_done()
            status = current_test.get_progress()
            version = current_test.get_version()
            return jsonify(running=True, percent=percent_done, status=status, version=version)    
            
    finally:
        schedule_lock.release()
    
    
@app.route('/finesse/get_<count>_logs', methods=['POST'])
def finesse_get_log(count):
    [out,err] = utils.git("--git-dir ./finesse_src/.git pull")
    
    [out,err] = utils.git("--git-dir ./finesse_src/.git log --max-count={0} --pretty=oneline".format(count))
    
    log_entries = out.split("\n")
    
    log2send = list()
    
    for e in log_entries:
        
        vals = e.split(" ",1)
        
        if len(vals[0]) > 8:
            if len(vals) > 1:
                message = vals[1]
            else:
                message = "[no commit message]"
            
            log2send.append({'commit':vals[0][0:8], 'message':message})
            
    
    return jsonify(logs=log2send)
    
class FinesseTestProcess(Thread):
    
    queue_time = None
    status = "Not started"
    built = False
    total_files = 100
    done_files = 0
    suite = ""
    git_commit = ""
    test_id = -1
    
    def __init__(self, git_commit, test_id, *args, **kqwargs):
        Thread.__init__(self)
        self.git_commit = git_commit
        self.queue_time = datetime.now()
        self.test_id = test_id
        
    def run(self):
        global current_test, scheduled_tests, schedule_lock
        
        for x in range(0,self.total_files):
            sleep(0.1)
            self.done_files += 1
            
            if x > 10:
                self.built = True
        
        # once done check if any other tests need to be ran
        schedule_lock.acquire()
        
        if len(scheduled_tests) > 0:
            current_test = scheduled_tests.pop(0)
            current_test.start()
        else:
            current_test = None
        
        schedule_lock.release()
        
        
    def percent_done(self):
        return 100.0*float(self.done_files)/float(self.total_files)
        
    def get_version(self):
        return self.git_commit
        
    def get_progress(self):
        if self.built:
            return '{0} out of {1} ({2})'.format(self.done_files, self.total_files, self.suite)
        else:
            return 'Building FINESSE...'
            
            
class FinesseTestProcessRun():
    def __init__(self, init_class=FinesseTestProcess):
        self.init_class = init_class    
        print "Init.."
        
    def start(self, *args, **kwargs):
        print "Calling..."
        fjp = self.init_class(*args, **kwargs)

        print '%s threaded process beginning.' % fjp.__class__.__name__
        print fjp.get_progress()

        fjp.start()

        while fjp.is_alive():
            
            sleep(1)
            
            print fjp.get_progress()

        print fjp.get_progress()
        print '%s threaded process complete. Now exiting.' % fjp.__class__.__name__
        
        
        
if __name__ == '__main__':
    
    # need local copy of src
    if not os.path.exists("./finesse_src"):
        print "finesse src folder didn't exist, cloning now..."
        utils.git("clone git://gitmaster.atlas.aei.uni-hannover.de/finesse/src.git finesse_src")
    else:
        # get the latest version for logs etc.
        utils.git("--git-dir ./finesse_src/.git pull")
        
    app.run(debug=True)