#!/usr/bin/env python

from pylab import *
from random import randint

groups = [	[ 'A', 'B', 'C', 'D', 'E' ],
			[ 'F', 'G', 'H', 'I' ],
			[ 'M', 'N', 'O', 'P', 'Q' ] ]

f = figure(1, figsize=(1,120))
for i in range(len(groups)):
	group = groups[i]
	subplot(len(groups), 1, i + 1)
	for treatment in group:
		plot(range(10), [ randint(1, 10) for x in range(10) ])

show()

