#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

N = 2

canexons = [ 9, 1 ]
altexons = [ 8, 0.5 ]
introns = [ 6, 0.2 ]

width = 0.5       # the width of the bars
leftspacing = 0.25 # space before first bar
padding = width * N+1

ind = (np.arange(N) * padding) + leftspacing # the x locations for the groups


fig = plt.figure()
plt.subplot(2,1,1)

rects2 = plt.bar(ind, altexons, width, color='y')
rects1 = plt.bar(ind+width, canexons, width, color='r')
rects3 = plt.bar(ind+2*width, introns, width, color='b')
plt.xticks(ind+1.5*width, ('ON', 'OFF'))
plt.legend( (rects1[0], rects2[0], rects3[0]), ('Canonical Exons', 'Alternative Exons', 'Introns') )
plt.title("First!")

plt.subplot(2,1,2)
rects2 = plt.bar(ind, altexons, width, color='y')
rects1 = plt.bar(ind+width, canexons, width, color='r')
rects3 = plt.bar(ind+2*width, introns, width, color='b')
plt.xticks(ind+1.5*width, ('ON', 'OFF'))
plt.legend( (rects1[0], rects2[0], rects3[0]), ('Canonical Exons', 'Alternative Exons', 'Introns') )
plt.title("Second!")


plt.show()

