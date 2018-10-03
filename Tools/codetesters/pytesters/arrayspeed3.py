#!/usr/bin/env python

from numpy import *
from pylab import rand
from time import clock
from scipy.stats.stats import nanmedian

import bp

def nan_median(a, axis=0):
	from numpy import ma
	masked_a = ma.masked_where(isnan(a), a)
	return ma.median(masked_a, axis=axis)

def my_median(vallist):
	num_vals = len(vallist)
	if num_vals == 0:
		return nan
	vallist.sort()
	if num_vals % 2 == 1: # odd
		index = (num_vals - 1) / 2
		return vallist[index]
	else: # even
		index = num_vals / 2
		return (vallist[index] + vallist[index - 1]) / 2

numtests = 100
testsize = 1000
pointlen = 3

t0 = clock()
natests = rand(numtests,testsize,pointlen)
# have to start with inf because list.remove(nan) doesn't remove nan
natests[natests > 0.9] = inf
tests = natests.tolist()
natests[natests==inf] = nan
for test in tests:
	for point in test:
		while inf in point:
			point.remove(inf)
t1 = clock()
print "list build time:", t1-t0

print natests

allmedians = []
t0 = clock()
for test in tests:
	medians = [ my_median(x) for x in test ]
	allmedians.append(medians)
t1 = clock()
print "list median time:", t1-t0

t0 = clock()
namedians = []
for natest in natests:
	thismed = []
	for point in natest:
		maskpoint = point[negative(isnan(point))]
		if len(maskpoint) > 0:
			med = median(maskpoint)
		else:
			med = nan
		thismed.append(med)
	namedians.append(thismed)
t1 = clock()
print "array nanmedian time:", t1-t0


t0 = clock()
namedians = []
for natest in natests:
	thismed = []
	for point in natest:
		if isnan(point).any():
			maskpoint = point[negative(isnan(point))]
			if len(maskpoint) > 0:
				med = median(maskpoint)
			else:
				med = nan
		else:
			med = median(point)
		thismed.append(med)
	namedians.append(thismed)
t1 = clock()
print "mixed median time:", t1-t0

t0 = clock()
namediansnew = []
for natest in natests:
	namedians.append(nan_median(natest, axis=1))
t1 = clock()
print "masked nanmedian time:", t1-t0

bp.breakpoint()

