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
from pykat.testing import test as finesse_test
import shutil
from pykat.testing.web import app

import os

global current_test, scheduled_tests, schedule_lock

test_id = 0
current_test = None
scheduled_tests = []
schedule_lock = Lock()
watcher = None

print "loading web interface"

class FinesseProcessWatcher(Thread):
    process_to_watch = None
    
    def __init__(self,*args, **kqwargs):
        Thread.__init__(self)
    
    def setProcessToWatch(self, process):
        print type(process)
        #if type(self.process_to_watch) is not finesse_test.FinesseTestProcess:
        #    raise Exception("Tried to watch something which wasn't a FinesseTestProcess")
    
        self.process_to_watch = process
        
    def run(self):
        global schedule_lock,current_test,scheduled_tests, watcher
        
        if self.process_to_watch is None:
            return
        
        #if type(self.process_to_watch) is not finesse_test.FinesseTestProcess:
        #    raise Exception("Tried to watch something which wasn't a FinesseTestProcess")
    
        print "Watcher is watching", self.process_to_watch
    
        self.process_to_watch.join()
        
        print "Watcher is continuing"
        
        # once done check if any other tests need to be ran
        schedule_lock.acquire()
        
        if len(scheduled_tests) > 0:
            print "Watcher starting next test"
            current_test = scheduled_tests.pop(0)
            watcher = FinesseProcessWatcher()
            watcher.setProcessToWatch(current_test)
            watcher.start()
            current_test.start()
        else:
            print "Watcher found no more tests to run"
            current_test = None
        
        schedule_lock.release()
        



@app.route('/')
def home_page():
    return render_template("finesse_test.html")

@app.route('/finesse/start_test', methods=["POST"])
def finesse_start_test():
    global current_test, test_id
    
    try:
        schedule_lock.acquire()
        
        test_id += 1
        
        git_commit = request.json['git_commit']
        
        print app.instance_path
        
        TEST_OUTPUT_PATH = os.path.join(app.instance_path, "tests")
        if not os.path.exists(TEST_OUTPUT_PATH):
            os.mkdir(TEST_OUTPUT_PATH)
            
        TEST_RUN_PATH = os.path.join(TEST_OUTPUT_PATH, str(test_id))
        print TEST_RUN_PATH
        
        if os.path.exists(TEST_RUN_PATH):
            shutil.rmtree(TEST_RUN_PATH)
            
        os.mkdir(TEST_RUN_PATH)
        
        test = finesse_test.FinesseTestProcess(os.path.join(app.instance_path, "finesse_test"), 
                                      TEST_RUN_PATH,
                                      git_commit, 
                                      run_fast=False, suites=[], test_id=test_id,
                                      emails="", nobuild=True)
        
        # check if anything is running and if it
        # isn't start this test off
        if current_test is None:
            print "running test"
            current_test = test
            watcher = FinesseProcessWatcher()
            watcher.setProcessToWatch(test)
            watcher.start()
            
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
                       