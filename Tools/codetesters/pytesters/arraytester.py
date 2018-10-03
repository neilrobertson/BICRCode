#!/usr/bin/env python

from pylab import rand
from numpy import median, nan, isnan, ma

numtests = 1
testsize = 100
pointlen = 3

test = rand(testsize,pointlen)

test[test > 0.9] = nan

print test
mtest = ma.masked_where(isnan(test), test)

print ma.median(mtest, axis=1)

