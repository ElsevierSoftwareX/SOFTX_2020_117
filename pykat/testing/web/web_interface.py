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
from pykat.testing import utils

from pykat.testing.web import app

import os

global current_test, scheduled_tests, schedule_lock

test_id = 0
current_test = None
scheduled_tests = []
schedule_lock = Lock()
    
print "loading web interface"

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
        
        