#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

N = 6

list1 = [ 1,2,1,5,2,4 ]
#list2 = [ 1,2,1,5,2,4 ]
list2 = [ 0.8,1.9,0.9,5,2.1,4.1 ]

width = 0.5       # the width of the bars
leftspacing = 0.25 # space before first bar
padding = width * 3

ind = (np.arange(N) * padding) + leftspacing # the x locations for the groups


fig = plt.figure()

print ind

rects1 = plt.bar(ind+width, list1, width, color='r')
rects2 = plt.bar(ind+2*width, list2, width, color='b')
plt.xticks(ind+2*width, ('e1', 'i1', 'e2', 'i2', 'e3', 'i3'))


plt.legend( (rects1[0], rects2[0]), ('Medians', 'Means') )

plt.show()

