from threading import Thread, Lock
from time import sleep
from uuid import uuid4
from flask import Flask
from flask import jsonify
from flask import Module
from flask import request, make_response
from flask import render_template, url_for
from datetime import datetime
from collections import namedtuple
from pykat.testing import utils
from pykat.testing import test as finesse_test
import shutil
from pykat.testing.web import app

from CodernityDB.database_thread_safe import ThreadSafeDatabase
from CodernityDB.database import RecordNotFound

from pykat.testing.web.database_indices import TestIDIndex, SrcCommitIndex

import os, sys

global current_test, scheduled_tests, schedule_lock

test_id = 0
current_test = None
scheduled_tests = []
schedule_lock = Lock()
watcher = None

print "Starting up database"
        
DB_PATH = os.path.join(app.instance_path,"db")

db = ThreadSafeDatabase(DB_PATH)

if db.exists():
    db.open()
    print db.get_db_details()
    print "Reindexing..."
    db.reindex()
    print "Done reindexing"
    
    # update the test_id from the previous tests that have been run
    for a in db.all('testid'):
        key = int(a['key'])
        if key > test_id:
            test_id = key
           
    print "Current test_id: " + str(test_id)
    
else:
    db.create()
    db.add_index(TestIDIndex(db.path, 'testid'))
    db.add_index(SrcCommitIndex(db.path, 'src_commit'))
    
    
print "loading web interface"

# should be called with the correct locks already
# applied
def __run_new(test):
    global current_test,watcher,scheduled_tests
    # check if anything is running and if it
    # isn't start this test off
    if current_test is None:
        print "running test"
        current_test = test
        # create watcher thread which will start the test
        watcher = FinesseProcessWatcher()
        watcher.setProcessToWatch(test)
        watcher.start()
    else:
        print "queuing test"
        scheduled_tests.append(test)

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
        start = datetime.now()
        
        self.process_to_watch.start()
        self.process_to_watch.join()
        
        print "Watcher is continuing"
        
        doc = db.get('testid', self.process_to_watch.test_id,with_doc=True)["doc"]
        doc["cancelled"] = self.process_to_watch.cancelling
        doc["error"] = self.process_to_watch.errorOccurred
        doc["startTime"] = str(start)
        doc["endTime"] = str(datetime.now())
        doc["testRun"] = self.process_to_watch.finished_test
        doc["diffFound"] = self.process_to_watch.diffFound
        
        db.update(doc)
        
        try:
            # once done check if any other tests need to be ran
            schedule_lock.acquire()
        
            if len(scheduled_tests) > 0:
                print "Watcher starting next test"
                current_test = scheduled_tests.pop(0)
                watcher = FinesseProcessWatcher()
                watcher.setProcessToWatch(current_test)
                watcher.start()
            else:
                print "Watcher found no more tests to run"
                current_test = None
                watcher = None
        finally:
            schedule_lock.release()
                

@app.route('/finesse/cancel_test_<id>', methods=["POST"])
def finesse_cancel_test(id):
    print id
    
    if int(id) >= 0:
        id = int(id)
        print "Cancelling " + str(id)
        
        try:
            # get lock here so that watcher doesn't interfere
            # with removing/starting new tests
            schedule_lock.acquire()
                        
            if current_test is not None:
                if current_test.test_id == id:
                    current_test.cancelling = True
                    print "Cancelling Current Test"
                    return str(id)
        
            ix = 0
            remove = -1;
            
            for t in scheduled_tests:
                if t.test_id == id:
                    remove = ix
                    break
                ix += 1
             
            if remove > -1:
                print "Cancelled queued test"
                scheduled_tests[remove].cancelling = True
                scheduled_tests.pop(remove)
                return str(id)
            
            print "Nothing cancelled"
            return "0"
        finally:
            schedule_lock.release()     

@app.route('/')
@app.route('/finesse/')
def home_page():
    return render_template("finesse_test.html")

@app.route('/finesse/rerun_test_<id>', methods=["POST"])
def finesse_start_rerun(id):
    id = int(id)
    try:
        schedule_lock.acquire()
        
        test_prev_doc = db.get('testid', id, with_doc=True)["doc"]
        
        if "rerun"  in test_prev_doc:
            test_prev_doc["rerun"] += 1
        else:
            test_prev_doc["rerun"] = 1
        
        db.update(test_prev_doc)
        
        TEST_OUTPUT_PATH = os.path.join(app.instance_path, "tests")
        
        if not os.path.exists(TEST_OUTPUT_PATH):
            raise Exception("Output path for rerun of test " + id + " does not exist")
        
        TEST_RUN_PATH_OLD = os.path.join(TEST_OUTPUT_PATH, str(id))        
        TEST_RUN_PATH = os.path.join(TEST_OUTPUT_PATH, str(id) + "_rerun_" + str(test_prev_doc["rerun"]))
        
        if not os.path.exists(TEST_RUN_PATH_OLD):
            NOBUILD = True
            shutil.copytree(os.path.join(TEST_RUN_PATH_OLD,"build"), os.path.join(TEST_RUN_PATH,"build"))
        else:
            NOBUILD = False
        
        test = finesse_test.FinesseTestProcess(os.path.join(app.instance_path, "finesse_test"), 
                                      TEST_RUN_PATH,
                                      test_prev_doc["git_commit"], 
                                      run_fast=True, suites=[], test_id=id,
                                      emails="", nobuild=NOBUILD)
        
        
        
        __run_new(test)
    finally:
        schedule_lock.release()

    return "ok"
    
@app.route('/finesse/start_test', methods=["POST"])
def finesse_start_test():
    global current_test, test_id
    
    try:
        schedule_lock.acquire()
        
        test_id += 1
        
        git_commit = request.json['git_commit']
                
        TEST_OUTPUT_PATH = os.path.join(app.instance_path, "tests")
        if not os.path.exists(TEST_OUTPUT_PATH):
            os.mkdir(TEST_OUTPUT_PATH)
            
        TEST_RUN_PATH = os.path.join(TEST_OUTPUT_PATH, str(test_id))
        
        if os.path.exists(TEST_RUN_PATH):
            shutil.rmtree(TEST_RUN_PATH)
            
        os.mkdir(TEST_RUN_PATH)
        
        test = finesse_test.FinesseTestProcess(os.path.join(app.instance_path, "finesse_test"), 
                                      TEST_RUN_PATH,
                                      git_commit, 
                                      run_fast=True, suites=[], test_id=test_id,
                                      emails="", nobuild=False)
        
        db.insert(dict(t="test",
                       test_id=test.test_id,
                       git_commit=test.get_version(),
                       cancelled=test.cancelling,
                       error=test.errorOccurred,
                       diffFound=test.diffFound,
                       testRun=test.finished_test))
                       
        __run_new(test)
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
            test_id = current_test.test_id
            cancelling = current_test.cancelling
            percent_done = current_test.percent_done()
            status = current_test.get_progress()
            version = current_test.get_version()
            return jsonify(cancelling=cancelling, running=True, id=test_id, percent=percent_done, status=status, version=version)    
            
    finally:
        schedule_lock.release()
    
@app.route('/finesse/get_branches', methods=["POST"])
def finesse_get_branches():
    os.chdir(os.path.join(app.instance_path,"finesse_src"))
        
    try:
        [out,err] = utils.git("branch -a")
    except Exception as ex:
        print "git branch error : " + str(ex)
    
    branches = list()
    
    for b in out.split("\n"):
        vals = b.split("/")
        if len(vals) >= 3:
            branches.append(vals[2].split(" ")[0])

    return jsonify(branches=branches)
    
@app.route('/finesse/get_<count>_<branch>_logs', methods=['POST'])
def finesse_get_log(count,branch):
    os.chdir(os.path.join(app.instance_path,"finesse_src"))
        
    print "!!!!", count, branch
    try:
        [out,err] = utils.git("checkout " + branch)
        [out,err] = utils.git("pull")
    except Exception as ex:
        print "git pull error : " + str(ex)
    
    [out,err] = utils.git("log --max-count={0} --pretty=oneline".format(count))
    
    log_entries = out.split("\n")
    
    log2send = list()
    
    for e in log_entries:
        
        vals = e.split(" ",1)
        vals[0] = vals[0].strip()
        
        if len(vals[0]) == 40:
            if len(vals) > 1:
                message = vals[1]
            else:
                message = "[no commit message]"
            
            log2send.append({'commit':vals[0], 'message':message})
            
    
    return jsonify(logs=log2send)
      
@app.route('/finesse/get_<count>_prev_tests', methods=["POST"])
def finesse_get_prev_tests(count):
    global test_id
    
    rtn = list()
    max = test_id
    min = max - int(count)

    if min < 0:
        min = 0
    if min > max:
        min = max
    
    try:
        data = db.all('testid',with_doc=True)
        
        for a in data:
            
            i = a["doc"]
            
            err = (not i['error'] is None)
                        
            if "startTime" in i:
                startTime = i["startTime"]
            else:
                startTime = ""
            
            if "endTime" in i:
                endTime = i["endTime"]
            else:
                endTime = ""
            
            global current_test
            
            if current_test is not None and current_test.test_id == i["test_id"]:
                status = "Running"
            elif err:
                status = "Test Exception"
            elif i["cancelled"] == True:
                status = "Cancelled"
            elif "diffFound" in i and i["diffFound"] == True:
                status = "ERRORS"
            elif "testRun" in i and i["testRun"] == False:
                status = "Not started"
            else:
                status = "OK"
            
            obj = dict(test_id=i['test_id'],
                           git_commit=i['git_commit'],
                           status=status,
                           startTime=startTime,
                           endTime=endTime)
            
            rtn.insert(0,obj)
           
        return jsonify(tests=rtn)
        
    except RecordNotFound:
        print "exception"
        return jsonify(test=rtn)

@app.route('/finesse/view/<view_test_id>/<log>.log')
def finesse_view_make(view_test_id, log):
    view_test_id = int(view_test_id)
    
    if view_test_id > 0:
        TEST_PATH = os.path.join(app.instance_path, "tests", str(view_test_id))
        response = ""
        
        log_string = ""
        
        if log == "build":
            log = os.path.join(TEST_PATH, "build", "build.log")
            log_string = open(log).read()
        elif log == "make":
            log = os.path.join(TEST_PATH, "build", "make.log")
            log_string = open(log).read()
            
        response = make_response(log_string)
        response.headers["Content-type"] = "text/plain"
            
        return response
        
    else:
        return ""
@app.route('/finesse/view/<view_test_id>/', methods=["GET"])
def finesse_view(view_test_id):
    
    view_test_id = view_test_id
    
    return render_template("finesse_test_view.html",
                            view_test_id = view_test_id)

    
    
    
    
        