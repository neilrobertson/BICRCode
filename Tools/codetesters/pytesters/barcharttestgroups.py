#!/usr/bin/env python

import numpy as np
import pylab as plt

from breakpoint import *

plt.rcParams['ytick.major.pad'] = 30

N = 2

canexons = [ 9, 1 ]
altexons = [ 8, 0.5 ]
introns = [ 6, 0.2 ]

cestds = [ 0.8, 0.1 ]
altstds = [ 0.7, 0.09 ]
instds = [ 0.5, 0.01 ]

width = 0.5       # the width of the bars
leftspacing = 0.25 # space before first bar
padding = width * N+1


assert(len(canexons) == len(altexons) == len(introns))
allplotdata = []
for i in range(len(canexons)):
	subplotdata = (canexons[i], altexons[i], introns[i])
	errordata = (cestds[i], altstds[i], instds[i])
	allplotdata.append((subplotdata, errordata))
	
print allplotdata

base = leftspacing


fig = plt.figure()
plt.subplots_adjust(wspace=0.7)

plotnum = 1
numplots = len(allplotdata)
for subplotdata in allplotdata:
	points, errors = subplotdata
	s = plt.subplot(1, numplots, plotnum)
	canexon, altexon, intron = points
	cestd, altstd, instd = errors
	rects2 = plt.bar(base, altexon, width, color='lightgreen', yerr=altstd)
	rects1 = plt.bar(base+width, canexon, width, color='yellow', yerr=cestd)
	rects3 = plt.bar(base+2*width, intron, width, color='lightblue', yerr=instd)
	

#	plt.xticks(base+1.5*width, ('ON', 'OFF'))




	plt.legend( (rects1[0], rects2[0], rects3[0]), ('Canonical Exons', 'Alternative Exons', 'Introns') )
	plotnum += 1
	s.axes.set_xticklabels([])
	s.axes.set_xticks([])

plt.show()

