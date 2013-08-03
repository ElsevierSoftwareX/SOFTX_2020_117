#!/bin/python
from threading import Thread, Lock
from time import sleep
from optparse import OptionParser
import os
import subprocess as sub
import numpy as np
import difflib
from StringIO import StringIO
import shutil
import smtplib
import string
import time
import pickle
from datetime import datetime
from pykat.testing import utils
import sys, traceback

class DiffException(Exception):
	def __init__(self, msg, outfile):
		self.msg = msg
		self.outfile = outfile

class FinesseTestProcess(Thread):
    
    queue_time = None
    status = ""
    built = False
    total_kats = 0
    done_kats = 0
    git_commit = ""
    test_id = -1
    finished_test = False
    diff_rel_eps = 1e-13
    running_kat = ""
    running_suite = ""
    cancelling = False
    errorOccurred = None
    diffFound = False
    diffing = False
    kats_to_run = []
    
    # Define storage structures for generating report later
    kat_run_exceptions = {}
    output_differences = {}
    run_times = {}
        
    def __init__(self, TEST_DIR, BASE_DIR, test_commit, 
                 run_fast=False, kats=[], test_id="0",
                 git_bin="",emails="", nobuild=True,*args, **kqwargs):
        
        Thread.__init__(self)
        self.git_commit = test_commit
        
        if test_commit is None:
            raise Exception("A git commit ID must be provided for the test")
        
        self.queue_time = datetime.now()
        self.test_id = test_id
        self.TEST_DIR = TEST_DIR
        self.BASE_DIR = BASE_DIR
        self.emails = ""
        
        if type(nobuild) is str:
        
            if nobuild.lower() == "true":
                self.nobuild = True
            elif nobuild.lower() == "false":
                self.nobuild = False
            else:
                raise Exception("nobuild is not a boolean value")
                
        elif type(nobuild) is bool:
            self.nobuild = nobuild
        else:
            raise Exception("nobuild is not a boolean value")
        
        if type(run_fast) is str:
            if run_fast.lower() == "true":
                self.run_fast = True
            elif run_fast.lower() == "false":
                self.run_fast = False
            else:
                raise Exception("run_fast is not a boolean value")
                
        elif type(run_fast) is bool:
            self.run_fast = run_fast
        else:
            raise Exception("nobuild is not a boolean value")
            
        if not os.path.isdir(self.BASE_DIR):
            raise Exception("BASE_DIR was not a valid directory")
        
        if not os.path.isdir(self.TEST_DIR):
            raise Exception("TEST_DIR was not a valid directory, should point to a clone of the FINESSE test repository")
            
        self.kats_to_run = kats
        
        self.GIT_BIN = git_bin
    
    def cancelCheck(self):
        if self.cancelling:
            raise SystemExit()
    
    def percent_done(self):
        if self.total_kats == 0:
            return 0.0
        else:
            return 100.0*float(self.done_kats)/float(self.total_kats)
        
    def get_version(self):
        return self.git_commit
        
    def get_progress(self):
        if self.diffing:
            return 'Diffing {0} out of {1} ({2} in {3})'.format(self.done_kats, self.total_kats/2, self.running_kat, self.running_suite)
        if self.built:
            return 'Running {0} out of {1} ({2} in {3})'.format(self.done_kats, self.total_kats/2, self.running_kat, self.running_suite)
        else:
            return 'Building FINESSE executable'
            
    def startFinesseTest(self):
        if sys.platform == "win32":
            EXE = ".exe"
        else:
            EXE = ""
            
        self.built = False

        BUILD_PATH = os.path.join(self.BASE_DIR, "build")
                    
        # Firstly we need to build the latest version of finesse
        if not self.nobuild:
            print "deleting build dir..." + BUILD_PATH
            if os.path.exists(BUILD_PATH):
                shutil.rmtree(BUILD_PATH)

            print "Checking out finesse base..."
            utils.git(["clone","git://gitmaster.atlas.aei.uni-hannover.de/finesse/base.git",BUILD_PATH])

            print "Checking out and building develop version of finesse " + self.git_commit
            
            SRC_PATH = os.path.join(BUILD_PATH,"src")
            
            if sys.platform == "win32":
                utils.runcmd(["bash","./finesse.sh","--checkout"],cwd=BUILD_PATH)
                self.cancelCheck()
                
                utils.git(["checkout",self.git_commit],cwd=SRC_PATH)
                self.cancelCheck()
                
                utils.runcmd(["bash","./finesse.sh","--build"],cwd=BUILD_PATH)
                self.cancelCheck()
            else:
                utils.runcmd(["./finesse.sh","--checkout","develop"],cwd=BUILD_PATH)
                self.cancelCheck()
                
                utils.git(["checkout",self.git_commit],cwd=SRC_PATH)
                self.cancelCheck()
                
                utils.runcmd(["./finesse.sh","--build"],cwd=BUILD_PATH)
                self.cancelCheck()
                            
        FINESSE_EXE = os.path.join(self.BASE_DIR,"build","kat" + EXE)
        
        # check if kat runs
        if not os.path.exists(FINESSE_EXE):
            raise Exception("Kat file was not found in " + FINESSE_EXE)
        
        self.built = True
        
        print "kat file found in " + FINESSE_EXE
        
        OUTPUTS_DIR = os.path.join(self.BASE_DIR,"outputs")
        
        if os.path.isdir(OUTPUTS_DIR):
            print "deleting outputs dir..."
            shutil.rmtree(OUTPUTS_DIR)
            
        os.mkdir(OUTPUTS_DIR)
        
        os.environ["KATINI"] = os.path.join(self.TEST_DIR,"kat_test","kat.ini")
        
        self.cancelCheck()
        # Clean up and pull latest test repository
        print "Cleaning test repository..."
        utils.git(["clean","-xdf"], cwd=self.TEST_DIR)
        self.cancelCheck()
        utils.git(["reset","--hard"], cwd=self.TEST_DIR)
        self.cancelCheck()
        print "Pulling latest test..."
        utils.git(["pull"],cwd=self.TEST_DIR)
        self.cancelCheck()
    
        self.total_kats = 0
        
        # create dictionary structures
        # and count up total number of files to process
        for suite in self.kats_to_run.keys():
            self.kat_run_exceptions[suite] = {}
            self.output_differences[suite] = {}
            self.run_times[suite] = {}
            
            self.total_kats += len(self.kats_to_run[suite])
        
        # multiply as we include the diffining in the percentage
        # done
        self.total_kats *= 2
        
        for suite in self.kats_to_run.keys():
            self.cancelCheck()
            print "Running suite: " + suite + "..."
            kats = self.kats_to_run[suite]
            SUITE_PATH = os.path.join(self.TEST_DIR,"kat_test",suite)

            SUITE_OUTPUT_DIR = os.path.join(OUTPUTS_DIR,suite)
            os.mkdir(SUITE_OUTPUT_DIR)

            self.running_suite = suite
            
            for kat in kats:
                self.cancelCheck()
                self.running_kat = kat
                
                print self.get_progress()
                basename = os.path.splitext(kat)[0]

                if self.run_fast and ('map ' in open(kat).read()):
                    print "skipping " + kat			
                else:
                    try:
                        start = time.time()
                        out,err = utils.runcmd([FINESSE_EXE, "--noheader", kat],cwd=SUITE_PATH)
                        
                        OUT_FILE = os.path.join(SUITE_PATH,basename + ".out")
                        shutil.move(OUT_FILE, SUITE_OUTPUT_DIR)
                    
                    except utils.RunException as e:
                    
                        print "STDERR: " + e.out
                        print "STDOUT: " + e.err
                        
                        print "Error running " + kat + ": " + e.err
                        self.kat_run_exceptions[suite][kat] = e
                        # if this happens a difference definitely is found
                        self.diffFound = True
                    finally:
                        self.run_times[suite][kat] = time.time()-start
                        self.done_kats += 1

        self.cancelCheck()
        
        for suite in self.kats_to_run.keys():
            if len(self.kat_run_exceptions[suite].keys()) > 0:
                print "Could not run the following kats:\n" + "\n".join(self.kat_run_exceptions[suite].keys()) + " in " + suite
            else:
                print "No errors whilst running" + suite

        self.diffing = True
        
        # Now we have generated the output files compare them to the references
        for suite in self.kats_to_run.keys():
            self.cancelCheck()
            print "Diffing suite: " + suite + "..."

            outs = []
            SUITE_PATH = os.path.join(OUTPUTS_DIR,suite)

            for files in os.listdir(SUITE_PATH):
                if files.endswith(".out"):
                    outs.append(files)

            REF_DIR = os.path.join(self.TEST_DIR,"kat_test",suite,"reference")

            if not os.path.exists(REF_DIR):
                raise Exception("Suite reference directory doesn't exist: " + REF_DIR)
                
            for out in outs:
                self.cancelCheck()
                self.running_kat = out
                
                ref_file = os.path.join(REF_DIR,out)
                out_file = os.path.join(SUITE_PATH,out)
                
                if not os.path.exists(ref_file):
                    raise DiffException("Reference file doesn't exist for " + out, out)
                    
                ref_arr = np.loadtxt(ref_file)
                out_arr = np.loadtxt(out_file)

                if ref_arr.shape != out_arr.shape:
                    raise DiffException("Reference and output are different shapes", out)

                # for computing relative errors we need to make sure we
                # have no zeros in the data
                ref_arr_c = np.where(ref_arr == 0, ref_arr, 1)
                ref_arr_c[ref_arr_c==0] = 1

                rel_diff = np.abs(out_arr-ref_arr)/np.abs(ref_arr_c)

                diff = np.any(rel_diff >= self.diff_rel_eps)
                
                if diff:
                    self.diffFound = True
                    
                    # store the rows which are different
                    ix = np.where(rel_diff >= self.diff_rel_eps)[0][0]
                
                    self.output_differences[suite][out] = (ref_arr[ix],
                                                           out_arr[ix],
                                                           np.max(rel_diff))
                else:
                    max = np.max(rel_diff) 
                    
                    self.output_differences[suite][out] = ([],
                                                           [],
                                                           max)
                self.done_kats += 1
                
        REPORT_PATH = os.path.join(self.BASE_DIR,"reports")
        
        if not os.path.exists(REPORT_PATH):
            os.mkdir(REPORT_PATH)

        today = datetime.utcnow()
        reportname = os.path.join(REPORT_PATH,"report.log") #today.strftime('%d%m%y')
        print "Writing report to " + reportname
        
        f = open(reportname,'w')
        f.write("Python Nightly Test\n")
        f.write(today.strftime('%A, %d. %B %Y %I:%M%p') + "\n")

        # add kat file header
        p = sub.Popen([FINESSE_EXE], stdout=sub.PIPE, stderr=sub.PIPE)
        out, err = p.communicate()
        f.write(out)
        
        # Now time to generate a report...
        np.set_printoptions(precision=16)
        
        isError = False

        for suite in self.kats_to_run.keys():
            f.write("\n\n" + str(len(self.output_differences[suite].keys())) + " differences in suite " + suite)
            for k in self.output_differences[suite].keys():
                isError = True
                f.write(k + ":\n")
                f.write("     ref: " + str(self.output_differences[suite][k][0]) + "\n")
                f.write("     out: " + str(self.output_differences[suite][k][1]) + "\n")
                f.write("     Max relative difference: " + str(self.output_differences[suite][k][2]) + "\n")

            f.write("\n\n" + str(len(self.output_differences[suite].keys())) + " errors in suite " + suite)
            for k in self.kat_run_exceptions[suite].keys():
                isError = True
                f.write(k + ":\n")
                f.write("err: " + self.kat_run_exceptions[suite][k].err + "\n")

        f.close()
        
        #if self.emails:
        #    
        #    if isError:
        #        subject = "Finesse test ERROR"
        #    else:
        #        subject = "Finesse test OK"
        #    emails = self.emails

        #    args = ["mailx", "-s", subject, emails]
        #    p = sub.Popen(args, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE)
        #    r = open(reportname,"r")
        #    out, err = p.communicate(r.read())
        #else:
        #    print "No emails specified"

            
            
    def run(self):
        
        try:
            self.startFinesseTest()
        except Exception as ex:
            
            exc_type, exc_value, exc_traceback = sys.exc_info()
            
            self.errorOccurred = dict(value=str(exc_value), traceback=str(traceback.format_exc(5)))
            
            print "*** Exception for test_id = " + str(self.test_id)
            traceback.print_exception(exc_type, exc_value, exc_traceback,
                                      limit=5, file=sys.stdout)
        finally:
            self.finished_test = True
        
        

if __name__ == "__main__":
    
    parser = OptionParser()
    
    parser.add_option("-t","--test-dir",dest="test_dir",help="")
    parser.add_option("-b","--base-dir",dest="base_dir",help="")
    parser.add_option("-c","--test-commit",dest="test_commit",help="")
    parser.add_option("-s","--suites",dest="suites",help="comma delimited list of each suite to run")
    parser.add_option("-g","--git-bin",dest="git_bin", default="git",help="")
    parser.add_option("-e","--emails",dest="emails", help="")
    parser.add_option("-n","--no-build",default="False",dest="nobuild",action="store_true")
    parser.add_option("-f","--fast",default="True",dest="fast",action="store_true")

    options, args = parser.parse_args()

    if options.test_dir is None:
        print "--test-dir argument is missing"
        exit()
        
    if options.test_commit is None:
        print "--test-commit argument is missing"
        exit()
    
    if options.base_dir is None:
        options.base_dir = os.getcwd()
    
    if options.suites is None:
        suites = []
    else:
        suites = options.suites.split(",")
        
    test = FinesseTestProcess(options.test_dir,
                              options.base_dir,
                              options.test_commit,
                              run_fast=options.fast,
                              suites=suites,
                              git_bin=options.git_bin,
                              emails=options.emails,
                              nobuild=options.nobuild)
    test.run()