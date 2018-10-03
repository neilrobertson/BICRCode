#!/usr/bin/env python

from pylab import *
x = [1,2,2,3,5,6,7,8,9]
y = [2,3,4,3,4,4,8,9,9]
z = [9,9,9,8,7,7,6,5,5]
figure(1)
plot(x,y,x,z)
legend(["y", "z"], loc="upper left")
#savefig("output.png")
show()


