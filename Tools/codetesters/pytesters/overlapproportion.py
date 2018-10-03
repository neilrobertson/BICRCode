#!/usr/bin/env python

# start end end is the bin
# lb and ub are the interval 
# we return how much of the interval overlaps the bin
def howMuchOverlap(start, end, lb, ub):
	binsize = end - start
	samplelength = ub - lb
	# outside - not in bin at all
	if ub < start or lb > end:
		return 0
	# completely inside bin
	if ub <= end and lb >= start:
		return 1
	# surrounding bin
	if ub >= end and lb <= start:
		return float(binsize) / float(samplelength)
	# to left of bin
	if lb < start and ub > start and ub < end:
		inbin = ub - start
		return float(inbin) / float(samplelength)
	# to right of bin
	if ub > end and lb < end and lb > start:
		inbin = end - lb
		return float(inbin) / float(samplelength)	


if __name__ == "__main__":
	# an overlap proportion should be the same regardless of the relative size of each bin
	# Each pair should be the same
	# no overlap
	print howMuchOverlap(10, 20, 25, 30)
	print howMuchOverlap(25, 30, 10, 20)
	# completely inside	
	print howMuchOverlap(10, 30, 25, 35)
	print howMuchOverlap(25, 30, 10, 40)
	print howMuchOverlap(0, 25, 10, 40)
	print howMuchOverlap(30, 50, 10, 40)
	# some overlaps. 
	print howMuchOverlap(10, 40, 25, 30)
	print howMuchOverlap(25, 30, 10, 40)
