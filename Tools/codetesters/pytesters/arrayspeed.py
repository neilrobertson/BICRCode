#!/usr/bin/env python

from random import random
import numpy
import time
import math
from scipy.stats.stats import nanmedian

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

threshold = 1e-5
def listsEqual(list1, list2):
	if len(list1) != len(list2):
		return False
	for i in range(len(list1)):
		diff = list1[i] - list2[i]
		if diff < 0:
			diff = -diff
		if diff > threshold:
			return False
	return True

def mean(vallist):
	return sum(vallist) / len(vallist)
	
def median(vallist):
	num_vals = len(vallist)
	vallist.sort()
	if num_vals % 2 == 1: # odd
		index = (num_vals - 1) / 2
		return vallist[index]
	else: # even
		index = num_vals / 2
		return (vallist[index] + vallist[index - 1]) / 2

# Finds the median of all values
def computeMedian(ma_map):
	all_vals = ma_map.values()
	return median(all_vals)

		
# Centre values around the median
def normaliseToMedian(mylist, median):
	return [ val / median for val in mylist ]
	

def log2all(mylist):
	return [ math.log(val) / math.log(2.0) for val in mylist if val > 0 ]


#####
# script starts here
#####		

numtests = 1000
testsize = 5000
testlb = 1.0
testub = 10.0

t0 = time.clock()
# build lists in advance so we can test on them separately.
tests = []
# number of tests
for i in range(numtests):
	testarray = [ (random() * (testub - testlb))  + testlb for i in range(testsize) ]
	tests.append(testarray)	
t1 = time.clock()
print "time to generate arrays:", (t1-t0)

# validation records
means = []
medians = []
normalised = []
logged = []

nameans = []
namedians = []
nanormalised = []
nalogged = []

t0 = time.clock()
#####
# regular lists test
#####
for test in tests:	
	# find mean
	means.append(mean(test))
	# find median
	medianval = median(test)
	medians.append(medianval)
	# normalise to median
	norm = normaliseToMedian(test, medianval)
	normalised.append(norm)
	# log2
	logged.append(log2all(norm))
t1 = time.clock()
print "time for list test:", (t1-t0)

t0 = time.clock()
# convert to numpy arrays
natests = []
for test in tests:
	natest = numpy.array(test)
	natests.append(natest)
t1 = time.clock()
print "time to convert to numpy arrays:", (t1-t0)

#####
# python array test
#####

t0 = time.clock()
#####
# numpy array test
#####
for test in natests:	
	# find mean
	nameans.append(test.mean())
	# find median
	medianval = nanmedian(test)
	namedians.append(medianval)
	# normalise to median
	norm = test / medianval
	nanormalised.append(norm)
	# log2
	nalogged.append(numpy.log2(norm))
t1 = time.clock()
print "time for numpy test:", (t1-t0)
	
#####
# validation
#####
try:
	assert(listsEqual(means, nameans))
except AssertionError:
	print "means do not match:"
	print "means:", means
	print "nameans:", nameans
try:
	assert(listsEqual(medians, namedians))
except AssertionError:
	print "medians do not match:"
	print "medians:", medians
	print "namedians:", namedians
assert(False not in (normalised == numpy.array(nanormalised)))
assert(False not in (logged == numpy.array(nalogged)))
	
