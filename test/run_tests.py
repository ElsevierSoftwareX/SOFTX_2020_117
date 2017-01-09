import pykat
import traceback
import os
import sys
import warnings

warnings.filterwarnings('error')

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    
errors = []

testdir = os.getcwd()



for path, folders, files in os.walk("./test_scripts"):
    
    for filename in files:
        if filename.endswith(".py"):
            filename = os.path.join(path, filename)
            
            with open(filename) as f:
                print("!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
                print("RUNNING: " + filename)
                try:
                    os.chdir(path)
                    code = compile(f.read(), filename, 'exec')
                    exec(code)
                except Exception as ex:
                    print(bcolors.FAIL)
                    (_type, value, tb) = sys.exc_info()
                    
                    print("EXCEPTION: "+ str(value))
                    
                    print(traceback.format_exc())
                    
                    errors.append(filename)
                    print(bcolors.ENDC)
                    sys.stdout.flush()
                    sys.stderr.flush()
                finally:
                    os.chdir(testdir)
                print("!------------------------------------------------------------------------------------------------")
                print("")
                print("")
    

if len(errors) > 0:
    print("\nFAILED !!!\n")
    print("The following files failed: ")
    for e in errors:
        print(" - " + e)
    
    sys.exit(1)
else:
    print("\nPASSED\n")
    sys.exit(0)