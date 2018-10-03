#!/usr/bin/env python

import pp
import random
import time
from string import lowercase

ncpus = 3

def timedCharDump(waittime, char):
	time.sleep(waittime)
	mycore = open("/proc/%i/stat" % os.getpid()).read().split()[38]
	print "I'm doing stuff!", mycore, char
	return [ char ]
	
job_server = pp.Server(ncpus, ppservers=())

jobdetails = [ (random.random(), letter) for letter in lowercase ]

jobs = [ job_server.submit(timedCharDump,(jinput1, jinput2), (), ("os", "time",)) for jinput1, jinput2 in jobdetails ]

for job in jobs:
	print job()[0]
