import csv
import os
import sys

from bedgenerators import *

from bp import breakpoint

# this is the simplest entry point into these functions
# pass in two bed files and it will return a generator to give you overlaps between the two
# you can optionally pass in a tolerance for the overlaps 
# and a list of the files that are sorted "left" or "ritht" or both.
# note that this is separated from the generator function itself - overlapBedFilesPlaneSweepGen()
# this is because otherwise, file missing exceptions won't be thrown until we call next(), which is a bit confusing
def overlapBedFilesPlaneSweep(leftbedfile, rightbedfile, tolerance=0, sortedfiles=[]):
	assert os.access(leftbedfile, os.R_OK), "missing file: %s" % leftbedfile
	assert os.access(rightbedfile, os.R_OK), "missing file: %s" % rightbedfile
	return overlapBedFilesPlaneSweepGen(leftbedfile, rightbedfile, tolerance, sortedfiles)

def overlapBedFilesPlaneSweepGen(leftbedfile, rightbedfile, tolerance=0, sortedfiles=[]):
	if "left" in sortedfiles:
		print >> sys.stderr, "using presorted file:", leftbedfile
		leftgen = bedGenerator(leftbedfile, tolerance=tolerance)
	else:
		leftgen = sortBedGenerator(leftbedfile, tolerance=tolerance)
	if "right" in sortedfiles:
		print >> sys.stderr, "using presorted file:", rightbedfile
		rightgen = bedGenerator(rightbedfile)
	else:
		rightgen = sortBedGenerator(rightbedfile)
	ogps = overlapGeneratorPlaneSweep(leftgen, rightgen)
	# To support a tolerance, we need to add a thin wrapper here to undo the tolerance modification
	for pair in ogps:
		lefti, righti = pair
		if tolerance > 0:
			lefti = (lefti[0], lefti[1]+tolerance, lefti[2]-tolerance, lefti[3])
		yield (lefti, righti)

# wrapper on overlapBedFilesPlaneSweep() to group hits by the intervals on the left
# this is useful when you want to look at all the hits for one interval at once before you act on any one hit
def overlapBedFilesPlaneSweepByLeft(leftbedfile, rightbedfile, tolerance=0, sortedfiles=[]):
	og = overlapBedFilesPlaneSweep(leftbedfile, rightbedfile, tolerance=tolerance, sortedfiles=sortedfiles)
	lastleft = None
	left, right = og.next()
	# pair of (leftinterval, set) where the set contains all the right intervals that overlap left
	current = (left, set())
	current[1].add(right)
	for left, right in og:
		if left != lastleft:
			yield current
			current = (left, set())
		current[1].add(right)
		lastleft = left
	yield current

# uses a plane sweep to determine the overlaps between two sets of genomic intervals
# they must come in as iterators
# sorted, with the same chromosome order
def overlapGeneratorPlaneSweep(left, right):
	# rightcache holds the points we've retrieved from the right set
	# that match the current left set chromosome, but we haven't finished filing yet
	rightcache = []
	# nextcache holds points from right that from a later chromosome 
	# to what we're currently looking at in left
	nextcache = []

	lastchrom = None
	for linterval in left:
		lchrom, lstart, lend, lvalue = linterval

		# slight complex mechanics to handle a change in chromosome in left
		if lchrom != lastchrom:
			# still some left over in right
			# make sure there aren't any points leftover in right
			# from the last chromosome from left
			if rightcache != []:
				righti = rightcache[-1]
				while righti[0] < lchrom:
					righti = right.next()
				nextcache = [ righti ]
			rightcache = []
			assert len(nextcache) <= 1
			# the next cache is from a chromosome *before* the left chrom
			# suck up all the points from before left's chrom 
			# - they will never overlap now
			if nextcache != [] and nextcache[0][0] < lchrom:
				righti = nextcache[0]
				while righti[0] < lchrom:
					righti = right.next()
				assert righti[0] >= lchrom
				nextcache = [ righti ]
			# if the point in nextcache matches the left chrom
			# nextcache is now current, so make it rightcach
			if nextcache != [] and nextcache[0][0] == lchrom:
				rightcache = nextcache[:]
				nextcache = []
			assert rightcache == [] or rightcache[0][0] == lchrom
			
		# if the last point in the rightcache overlaps this left point,
		# there might be more points in the right that are relevant now
		while nextcache == [] and (rightcache == [] or rightcache[-1][1] <= lend):
			try:
				righti = right.next()
			except StopIteration:
				break
			if righti[0] == lchrom:
				rightcache.append(righti)
			elif righti[0] == lastchrom or (righti[0] == lchrom and righti[2] <= lstart):
				pass
			else:
				nextcache.append(righti)
			
		removecount = 0
		# the bit where we actually check and report the overlaps
		for rinterval in rightcache:
			# we've gone past this one. Remove it
			if rinterval[2] <= lstart:
				removecount += 1
			# we haven't reached this one yet. 
			# We must have filed everything relevant from right.
			elif rinterval[1] >= lend:
				break
			else:
				# if we get here, this is an overlap. Yield it.
				yield((linterval, rinterval))
		
		rightcache = rightcache[removecount:]
					
		lastchrom = lchrom

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: overlapgenerator.py <leftfile> <rightfile>"
		sys.exit(1)
	leftfile = sys.argv[1]
	rightfile = sys.argv[2]
# 	leftfile = "/home/pzs/sreenivas/restrictionenzymes/csp6i.mm9.bed.sort"
# 	rightfile = "testinput/random1000.bed.sort"
	og = overlapBedFilesPlaneSweep(leftfile, rightfile)
#	overlaps = [ overlap for overlap in og ] 
#	print overlaps
#	print len(overlaps)
	outfile = "ogtester2.txt"
	ofh = open(outfile, "w")
	for overlap in og:
		print >> ofh, overlap
