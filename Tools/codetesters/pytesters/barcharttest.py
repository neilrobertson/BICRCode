#!/usr/bin/env python

import numpy as np
import pylab as plt
from matplotlib.font_manager import FontProperties

from breakpoint import *

def set_ticks_both(axis):
	ticks = list( axis.majorTicks ) # a copy
	ticks.extend( axis.minorTicks )

	for t in ticks:
		t.tick1On = True # tick marker on left (or bottom)
		t.tick2On = True # tick marker on right (or top)
		t.label1On = True # tick label marker on left (or bottom)
		t.label2On = True # tick label on right (or top)

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

ind = (np.arange(N) * padding) + leftspacing # the x locations for the groups


fig = plt.figure()

rects2 = plt.bar(ind, altexons, width, color=(0.0,0.3,0.0), yerr=altstds)
rects1 = plt.bar(ind+width, canexons, width, color=(0.0,0.6,0.0), yerr=cestds)
rects3 = plt.bar(ind+2*width, introns, width, color=(0.0,0.9,0.0), yerr=instds)

plt.xticks(ind+1.5*width, ('ON', 'OFF'), size=5)

set_ticks_both(fig.gca().yaxis)


plt.legend( (rects1[0], rects2[0], rects3[0]), ('Canonical Exons', 'Alternative Exons', 'Introns') )

plt.show()

