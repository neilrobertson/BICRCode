#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

N = 3

canexons = [ 9, 1, 2 ]
altexons = [ 8, 0.5, 3 ]
introns = [ 6, 0.2 ]

cestds = [ 0.8, 0.1 ]
altstds = [ 0.7, 0.09 ]
instds = [ 0.5, 0.01 ]

width = 0.5       # the width of the bars
leftspacing = 0.25 # space before first bar
padding = width * N+1

ind = (np.arange(N) * padding) + leftspacing # the x locations for the groups


fig = plt.figure()


for value in altexons:
	rects2 = plt.bar(ind, value, width, color='y', yerr=altstds)
for value in canexons:
	rects1 = plt.bar(ind+width, value, width, color='r', yerr=cestds)
#rects3 = plt.bar(ind+2*width, introns, width, color='b', yerr=instds)
plt.xticks(ind+1.5*width, ('ON', 'OFF'))


plt.legend( (rects1[0], rects2[0]), ('Canonical Exons', 'Alternative Exons') )

plt.show()

