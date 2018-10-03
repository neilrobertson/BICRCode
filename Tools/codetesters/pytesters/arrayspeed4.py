#!/usr/bin/env python

from pylab import rand
from numpy import median, isnan
from time import clock

def my_median(vallist):
	if isnan(vallist).any():
		print "attempt to median list containing nan!"
		sys.exit(1)
	num_vals = len(vallist)
	vallist.sort()
	if num_vals % 2 == 1: # odd
		index = (num_vals - 1) / 2
		return vallist[index]
	else: # even
		index = num_vals / 2
		return (vallist[index] + vallist[index - 1]) / 2

samples = rand(100000, 3)


t0 = clock()
medians = [ my_median(sample) for sample in samples ]
t1 = clock()
print "list time:", t1-t0

t0 = clock()
namedians = [ median(sample) for sample in samples ]
t1 = clock()
print "numpy time:", t1-t0
