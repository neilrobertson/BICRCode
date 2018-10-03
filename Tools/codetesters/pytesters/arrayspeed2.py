#!/usr/bin/env python

from numpy import median, array
from pylab import rand
from time import clock

from resource import getrusage

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

numtests = 1000
numpoints = 5000
pointdepth = 3

def my_median(vallist):
	num_vals = len(vallist)
	vallist.sort()
	if num_vals % 2 == 1: # odd
		index = (num_vals - 1) / 2
		return vallist[index]
	else: # even
		index = num_vals / 2
		return (vallist[index] + vallist[index - 1]) / 2

t0 = clock()
natests = []
tests = []
# generate lists
for i in range(numtests):
	natest = rand(numpoints, pointdepth)
	test = natest.tolist()[:]
	natests.append(natest)
	tests.append(test)



allmedians = []

# list test
t0 = clock()
for test in tests:
	medians = [ my_median(vals) for vals in test ]
	allmedians.append(medians)
t1 = clock()
print "list time:", t1-t0


allnamedians = []

# numpy test
t0 = clock()
for natest in natests:
	namedians = median(natest.T)	
	allnamedians.append(namedians)
t1 = clock()
print "numpy time:", t1-t0


# validation
try:
	assert(False not in (array(allmedians) == array(allnamedians)))
except AssertionError:
	print allmedians
	print allnamedians

print getrusage(RUSAGE_SELF)
