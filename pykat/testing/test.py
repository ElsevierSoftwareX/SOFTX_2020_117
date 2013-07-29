 #!/bin/python

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
import datetime

options = None
diff_rel_eps = 1e-13
GIT_BIN = ""

class RunException(Exception):
	def __init__(self, returncode, args, err, out):
		self.returncode = returncode
		self.args = args
		self.err = err
		self.out = out

class DiffException(Exception):
	def __init__(self, msg, outfile):
		self.msg = msg
		self.outfile = outfile

def runcmd(args):
	p = sub.Popen(args, stdout=sub.PIPE, stderr=sub.PIPE)
	out, err = p.communicate()

	if p.returncode != 0:
		print err
		raise RunException(p.returncode, args, err, out)

	return [out,err]

def git(args):
	p = sub.Popen([GIT_BIN] + args, stdout=sub.PIPE, stderr=sub.PIPE)
	out, err = p.communicate()

	if p.returncode != 0:
		print err
		raise RunException(p.returncode, args, err, out)

	return [out, err]


	BASE_DIR = os.getcwd()

	if not options.suites:
		suites = ["physics","random"]				
	else:
		suites = []
		suites.extend(options.suites.split(","))


	if not options.git:
		GIT_BIN = "/usr/git/bin"
	else:
		GIT_BIN = options.git

	if not options.fast:
		run_fast = False
	else:
		run_fast = True
		print "Running fast test"

	# Firstly we need to build the latest version of finesse
	if os.path.isdir("build") and not options.nobuild:
		print "deleting build dir..."
		shutil.rmtree("build")

		print "Checking out finesse base..."
		git(["clone","git://gitmaster.atlas.aei.uni-hannover.de/finesse/base.git","build"])

		os.chdir("build")
		print "Checking out develop version of finesse..."
		runcmd(["./finesse.sh","--checkout","develop"])
		print "Building finesse..."
		runcmd(["./finesse.sh","--build"])
		os.chdir(BASE_DIR)

	# check if kat runs
	if not os.path.exists("./build/kat"):
		raise Exception("Kat file was not found")

	FINESSE_EXE = os.path.join(os.getcwd(),"build","kat")

	print "kat file found in " + FINESSE_EXE


	OUTPUTS_DIR = os.path.join(BASE_DIR,"outputs")
	if os.path.isdir(OUTPUTS_DIR):
		print "deleting outputs dir..."
		shutil.rmtree(OUTPUTS_DIR)

	os.mkdir(OUTPUTS_DIR)
	
	os.environ["KATINI"]=os.path.join(BASE_DIR,"build","kat.ini")
	
	# Clean up and pull latest test repository
	os.chdir(os.path.join(options.test_git))
	print "Cleaning test repository...."
	git(["clean","-xdf"])
	git(["reset","--hard"])
	print "Pulling latest test..."
	git(["pull"])

	# Define storage structures for generating report later

	kat_run_exceptions = {}
	output_differences = {}
	run_times = {}

	# create dictionary structures
	for suite in suites:
		kat_run_exceptions[suite] = {}
		output_differences[suite] = {}
		run_times[suite] = {}


	for suite in suites:
		print "Running suite: " + suite + "..."
		kats = []
		os.chdir(os.path.join(options.test_git,"kat_test",suite))

		for files in os.listdir("."):
			if files.endswith(".kat"):
				kats.append(files)

		SUITE_OUTPUT_DIR = os.path.join(OUTPUTS_DIR,suite)
		os.mkdir(SUITE_OUTPUT_DIR)

		for kat in kats:
			print "Running kat: " + kat
			basename = os.path.splitext(kat)[0]

			if run_fast and ('map ' in open(kat).read()):
				print "skipping " + kat			
			else:
				try:
					start = time.time()
					out,err = runcmd([FINESSE_EXE, "--noheader", kat])
					finish = time.time()-start
					run_times[suite][kat] = finish
					shutil.move(basename + ".out", SUITE_OUTPUT_DIR)
				except RunException as e:
					print "Error running " + kat + ": " + e.err
					kat_run_exceptions[suite][kat] = e

	for suite in suites:
		if len(kat_run_exceptions[suite].keys()) > 0:
			print "Could not run the following kats:\n" + "\n".join(kat_run_exceptions.keys()) + " in " + suite
		else:
			print "No errors whilst running" + suite

	
	# Now we have generated the output files compare them to the references
	for suite in suites:
		print "Diffing suite: " + suite + "..."

		outs = []
		os.chdir(os.path.join(OUTPUTS_DIR,suite))

		for files in os.listdir("."):
			if files.endswith(".out"):
				outs.append(files)

		REF_DIR = os.path.join(options.test_git,"kat_test",suite,"reference")

		if not os.path.exists(REF_DIR):
			raise Exception("Suite reference directory doesn't exist: " + REF_DIR)
		for out in outs:
			#print "Diffing " + out
			ref_file = os.path.join(REF_DIR,out)
			
			if not os.path.exists(ref_file):
				raise DiffException("Reference file doesn't exist for " + out, out)
			ref_arr = np.loadtxt(ref_file)
			out_arr = np.loadtxt(out)

			if ref_arr.shape != out_arr.shape:
				raise DiffException("Reference and output are different shapes", out)

			# for computing relative errors we need to make sure we
			# have no zeros in the data
			ref_arr_c = np.where(ref_arr == 0, ref_arr, 1)
			ref_arr_c[ref_arr_c==0] = 1

			rel_diff = np.abs(out_arr-ref_arr)/np.abs(ref_arr_c)

			diff = np.any(rel_diff >= diff_rel_eps)
			
			if diff:
				# store the rows which are different
				ix = np.where(rel_diff >= diff_rel_eps)[0][0]
				output_differences[suite][out] = (ref_arr[ix], out_arr[ix], np.max(rel_diff))


	os.chdir(BASE_DIR)
	if not os.path.exists("reports"):
		os.mkdir("reports")

	os.chdir("reports")
	today = datetime.datetime.utcnow()
	reportname = today.strftime('%d%m%y')
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

	for suite in suites:
		f.write("\n\n" + str(len(output_differences[suite].keys())) + " differences in suite " + suite)
		for k in output_differences[suite].keys():
			isError = True
			f.write(k + ":\n")
			f.write("     ref: " + str(output_differences[suite][k][0]) + "\n")
			f.write("     out: " + str(output_differences[suite][k][1]) + "\n")
			f.write("     Max relative difference: " + str(output_differences[suite][k][2]) + "\n")

		f.write("\n\n" + str(len(output_differences[suite].keys())) + " errors in suite " + suite)
		for k in kat_run_exceptions[suite].keys():
			isError = True
			f.write(k + ":\n")
			f.write("err: " + kat_run_exceptions[suite][k].err + "\n")
			
	

	f.close()
	
	if options.emails:
		
		if isError:
			subject = "Finesse test ERROR"
		else:
			subject = "Finesse test OK"

		emails = options.emails

		args = ["mailx", "-s", subject, emails]
		p = sub.Popen(args, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE)
		r = open(reportname,"r")
		out, err = p.communicate(r.read())
	else:
		print "No emails specified"


