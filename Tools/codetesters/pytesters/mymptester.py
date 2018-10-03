#!/usr/bin/env python

from multiprocessing import Pool, current_process
import random
import time
from string import lowercase

ncpus = 3

def timedCharDump(args):
	waittime, char = args
	time.sleep(waittime)
	
	print "I'm doing stuff!", current_process().name, char
	return char 
	
def timedCharDump2(waittime, char):
	time.sleep(waittime)
	
	print "I'm doing stuff!", current_process().name, char
	return char

pool = Pool(processes=ncpus)

jobdetails = [ (random.random(), letter) for letter in lowercase ]

times = [ j[0] for j in jobdetails ]
letters = [ j[1] for j in jobdetails ]

print pool.map(timedCharDump, jobdetails)

print pool.map(timedCharDump2, times, letters)
