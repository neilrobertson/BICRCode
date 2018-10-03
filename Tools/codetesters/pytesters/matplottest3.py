#!/usr/bin/env python

import pickle
from pylab import *

from matplotlib.font_manager import FontProperties

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

points = pickle.load(open("points.pic", "rb"))
positions = pickle.load(open("positions.pic", "rb"))

groups = [	["H2B", "H3", "histones" ], 
			[ "K4Me1", "K4Me3", "K27Me1", "K27Me3" ],
			[ "K9-K14ac", "K36Me1", "K36Me3", "K79Me1", "K79Me3" ]	]

xmin = min(positions)
xmax = max(positions)

axes = []

# compute axis bounds
for group in groups:
	ymax = max([ max(points[treatment]) for treatment in group ])
	ymin = min([ min(points[treatment]) for treatment in group ])
	axes.append([ int(100 * xmin) - 5, int(100 * xmax) + 5, ymin - 1, ymax + 6 ])


f = figure(1, figsize=(5,10))
subplots_adjust(hspace=0.7)
for i in range(len(groups)):
	group = groups[i]
	s = subplot(len(groups), 1, i+1)
	for treatment in group:
		try:
			plot([ int(x * 100) for x in positions ], points[treatment])
		except:
			pass
			breakpoint()
	s.axes.set_xticklabels([])
	s.axes.set_yticklabels([])
#	a = axis(axes[i])
#	a = s.add_axes(axes[i], xticklabels="")
#	a.xaxis.set_visible(False)
#	a.set_xticklabels([], visible=False)
# 	legend(group, loc="upper left", prop=FontProperties(size="smaller"))
# 	xlabel('Proportion along gene (%)')
# 	ylabel('Enrichment (zscored)')

#savefig("output.png")
show()
close(f)

