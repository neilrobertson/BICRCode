#!/usr/bin/env python

from numpy import linspace

mainsize = 9921

mainlist = range(mainsize)

numbreaks = 12

breakpositions = linspace(0, mainsize, numbreaks+1)

partitions = []
percent = 100 / numbreaks
for i in range(len(breakpositions)-1):
	thispos = int(breakpositions[i])
	nextpos = int(breakpositions[i+1])
	if thispos == 0:
		fromprop = 0
	else:
		fromprop = (thispos / float(mainsize)) * 100
	toprop = (nextpos / float(mainsize)) * 100
	print "%s, %s" % (round(fromprop,0), round(toprop,0))
	partitions.append(mainlist[thispos:nextpos])
	


print "number of partitions", len(partitions)
partlengths = [ len(x) for x in partitions ]
print "length of each partition", partlengths
print "total length", sum(partlengths)

assert len(partitions) == numbreaks
assert sum(partlengths) == mainsize
