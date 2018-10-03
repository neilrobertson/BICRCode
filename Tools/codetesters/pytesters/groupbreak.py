#!/usr/bin/env python

longlist = range(10)

breakpositions = range(0, 10, 4)
breaks = []
start = 0
for i in range(len(breakpositions) - 1):
	start, end = breakpositions[i], breakpositions[i+1]
	breaks.append((start, end))
breaks.append((breakpositions[-1], 10))

print breaks

brokenlist = [ longlist[r[0]:r[1]] for r in breaks ]

print brokenlist
