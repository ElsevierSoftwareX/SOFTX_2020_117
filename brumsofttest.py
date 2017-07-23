# this will work with python 2.7 ONLY!
import sys,urllib,json,time
from requests import get

ip = get('https://api.ipify.org').text
print 'My public IP address is:', ip

SERVER = sys.argv[1]
USER_ID = sys.argv[2]
TEST_PROFILE = sys.argv[3]
GIT_CHECKOUT = sys.argv[4]

submit_url = SERVER + 'cli/submit_test/user/' +\
			 USER_ID + '/testProfile/' + TEST_PROFILE +\
			 '/gitcheckout/' + GIT_CHECKOUT
print('Requesting Test on URL: ' + str(submit_url))
response = urllib.urlopen(submit_url)

if response.getcode() != 200:
	print(response.read())
	raise Exception()

status_url = SERVER + 'cli/get_status/session/' +\
			 response.read() + '/user/' + USER_ID

print('Getting status from url: ' + str(status_url))

while True:
	time.sleep(2)
	new_response = urllib.urlopen(status_url)
	status = json.load(new_response)['status']
	if status == "Failed":
		raise Exception("Test Failed!!")
	elif status == "Passed":
		exit()
