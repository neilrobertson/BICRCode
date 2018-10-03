#!/usr/bin/env python

import pylab

from bp import breakpoint

x = [ 5, 7, 11, 12, 13, 21, 31, 52, 62, 65, 73, 82, 91, 92, 93, 94 ]
x2 = [ 5, 7, 11, 12, 13, 21, 31, 52, 62, 65, 73, 82, 91, 92, 93, 94, 95 ]


n, bins, patches = pylab.hist([ x, x2 ], bins=(0,20,100), range=(0,100))

pylab.legend((patches[0][0], patches[1][0]), ("x", "x2"))
pylab.show()
