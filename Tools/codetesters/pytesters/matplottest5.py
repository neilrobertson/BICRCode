#!/usr/bin/env python
import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

from random import randint
from pylab import *

def set_ticks_both(axis):
	ticks = list( axis.majorTicks ) # a copy
	ticks.extend( axis.minorTicks )

	for t in ticks:
		t.tick1On = True # tick marker on left (or bottom)
		t.tick2On = True # tick marker on right (or top)
		t.label1On = True # tick label marker on left (or bottom)
		t.label2On = True # tick label on right (or top)

groups = [	[ 'A', 'B', 'C', 'D', 'E' ],
			[ 'F', 'G', 'H', 'I' ],
			[ 'M', 'N', 'O', 'P', 'Q' ] ]

x = range(9)
y = [ randint(0,10) for x0 in range(9) ]
z = [ randint(0,10) for x0 in range(9) ]
a = [ randint(0,10) for x0 in range(9) ]
b = [ randint(0,10) for x0 in range(9) ]

data = []
for group in groups:
	innerdata = []
	for treatment in group:
		thisdata = [ randint(0,10) for x0 in range(9) ]
		innerdata.append(thisdata)
	data.append(innerdata)
print data
		
f = figure(1)
for i in range(len(groups)):
	print "subplot(%s, 1, %s)" % (len(groups), i)
	subplot(len(groups), 1, i + 1)
	group = groups[i]
	for treatment in group:
		print treatment
		plot(x,y)
		plot(x,a)
		plot(x,b)

legend(("y", "a", "b"), loc="upper left")
xlabel('Proportion along gene')
ylabel('log2 chip enrichment')

s = subplot(212)
plot(x,z)
legend(("z"), loc="upper left")
xlabel('Proportion along gene')
ylabel('log2 chip enrichment')

set_ticks_both(s.yaxis)
#savefig("output.png")


show()
