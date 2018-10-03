#!/usr/bin/env python

# also validates the data
def clipRanges(regions):
	for i in range(len(regions) - 1):
		thispoint = regions[i]
		nextpoint = regions[i+1]
		assert thispoint[1] > thispoint[0] and	nextpoint[1] > nextpoint[0], "point read not valid"
		thisend = thispoint[1] 
		nextstart = nextpoint[0]
		diff = thisend - nextstart
		# a difference of zero is too close together
		if diff > -1:
			if diff % 2 == 1:
				diff += 1
			correction = diff / 2
			newend = thisend - correction
			newstart = newend + 1
			assert newend > thispoint[0] and nextpoint[1] > newstart, "new range not valid!"
			newthispoint = (thispoint[0], newend, thispoint[2])
			newnextpoint = (newstart, nextpoint[1], nextpoint[2])
			regions[i] = newthispoint
			regions[i+1] = newnextpoint
	return regions

regions = [ (0,10,2.5), (12,22,3.5), (15,25,1.2), (23, 30,0.01), (27, 37,1.23), (30, 35, 1.45) ]
regions2 = [ (0,10,2.5), (1,11,1.1), (2,12,1.2) ]

# works fine, produces [(0, 10, 2.5), (12, 18, 3.5), (19, 24, 1.2), (25, 28, 0.01), (29, 33, 1.23), (34, 35, 1.45)]
print clipRanges(regions)
# violates "new range not valid" assertion
print clipRanges(regions2)
