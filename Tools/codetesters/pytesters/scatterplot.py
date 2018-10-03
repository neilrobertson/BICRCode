#!/usr/bin/env python

import pylab 

f = pylab.figure(figsize=(40, 8))
for i in range(20):
	s = pylab.subplot(2, 10, i + 1)
	pylab.scatter([ 1, 2, 3], [1,2,3])
	pylab.title("plot " + str(i))
pylab.show()
