#!/usr/bin/env python

import numpy
from pylab import *

from breakpoint import *

x = [1,2,2,3,5,6,7,8,9]
y = [2,3,4,3,4,4,8,9,9]
z = [9,9,9,8,7,7,6,5,5]
a = [9,9,9,8,7,7,6,5,5]
b = [9,9,9,8,7,7,6,5,5]

f = figure(1)

p = plot(x,z, 'green')
#plot(x,z, 'green')
plot(x,y, 'darkblue')
l = legend(('z','y'))

#savefig("output.png")

f.axes[0].set_xticklabels([])
#breakpoint()
show()
