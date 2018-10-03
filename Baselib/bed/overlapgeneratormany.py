import csv
import os
import sys

from glob import glob

from bedgenerators import *

from bp import breakpoint

def overlapBedFilesPlaneSweepMany(leftbedfile, rightbedfilemap, tolerance=0, sortedfiles=[]):
	assert os.access(leftbedfile, os.R_OK), "missing file: %s" % leftbedfile
	for rightname in rightbedfilemap:
		rightbedfile = rightbedfilemap[rightname]
		assert os.access(rightbedfile, os.R_OK), "missing file: %s" % rightbedfile
	return overlapBedFilesPlaneSweepManyGen(leftbedfile, rightbedfilemap, tolerance, sortedfiles)

def overlapBedFilesPlaneSweepManyGen(leftbedfile, rightbedfilemap, tolerance=0, sortedfiles=[]):
	leftname = os.path.basename(leftbedfile)
	if leftname in sortedfiles:
		print >> sys.stderr, "using presorted file:", leftbedfile
		leftgen = bedGenerator(leftbedfile, tolerance=tolerance)
	else:
		leftgen = sortBedGenerator(leftbedfile, tolerance=tolerance)
	rightgens = {}
	for rightname in rightbedfilemap:
		rightbedfile = rightbedfilemap[rightname]
		if rightname in sortedfiles:
			print >> sys.stderr, "using presorted file:", rightbedfile
			rightgen = bedGenerator(rightbedfile)
		else:
			rightgen = sortBedGenerator(rightbedfile)
		rightgens[rightname] = rightgen
	ogpsm = overlapGeneratorPlaneSweepMany(leftgen, rightgens)
	# To support a tolerance, we need to add a thin wrapper here to undo the tolerance modification
	for pair in ogpsm:
		lefti, righti = pair
		if tolerance > 0:
			lefti = (lefti[0], lefti[1]+tolerance, lefti[2]-tolerance, lefti[3])
		yield (lefti, righti)

# this is just moving the code out of the main routine. 
# It works entirely in side-effects on the lists that are passed in.
def chromosomeChangeOver(lchrom, right, nextcache, rightcache):
	if rightcache != []:
		righti = rightcache[-1]
		# still some left over in right
		# make sure there aren't any points leftover in right
		# from the last chromosome from left
		startchrom = righti[0]
		while righti[0] < lchrom:
			righti = right.next()
			assert righti[0] >= startchrom, "unsorted file: %s" % (right)
		# the new chromosome must be further on than the last
		# otherwise the incoming file is not properly sorted
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
		assert righti[0] >= lchrom, "unsorted file: %s" % (right)
		nextcache = [ righti ]
	# if the point in nextcache matches the left chrom
	# nextcache is now current, so make it rightcache
	if nextcache != [] and nextcache[0][0] == lchrom:
		rightcache = nextcache[:]
		nextcache = []
	assert rightcache == [] or rightcache[0][0] == lchrom
	return rightcache, nextcache

def updatePoints(linterval, lastchrom, right, nextcache, rightcache):
	lchrom, lstart, lend, lvalue = linterval
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
			assert rightcache == [] or righti[0] > rightcache[0][0], "unsorted file: %s" % (right)
			nextcache.append(righti)
	return rightcache, nextcache


# uses a plane sweep to determine the overlaps between a base set of genomic intervals
# and a list of others
# they must come in as iterators
# sorted, with the same chromosome order
# rights is a map from a name for each right track and the right interval generator
def overlapGeneratorPlaneSweepMany(left, rights):
	# rightcaches holds the points we've retrieved from the right generators 
	# (dictionary, because there are multiple)
	# that match the current left set chromosome, but we haven't finished filing yet
	rightcaches = dict([ (rightname, []) for rightname in rights ])
#	rightcache = []
	# nextcache holds points from right files that from a later chromosome 
	# to what we're currently looking at in left
	nextcaches = dict([ (rightname, []) for rightname in rights ])
#	nextcache = []

	lastchrom = None
	for linterval in left:
		lchrom, lstart, lend, lvalue = linterval

		# slight complex mechanics to handle a change in chromosome in left
		for rightname in rights:
			rightcache = rightcaches[rightname]
			nextcache = nextcaches[rightname]
			right = rights[rightname]
			if lchrom != lastchrom:
				assert lastchrom < lchrom
				print "working at", lastchrom
				rightcache, nextcache = chromosomeChangeOver(lchrom, right, nextcache, rightcache)

			rightcache, nextcache = updatePoints(linterval, lastchrom, right, nextcache, rightcache)
			
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
					yield((linterval, (rightname,) + rinterval))

			rightcache = rightcache[removecount:]

			rightcaches[rightname] = rightcache
			nextcaches[rightname] = nextcache

		lastchrom = lchrom

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "usage: overlapgenerator.py <leftfile> <rightdir> <extension>"
		sys.exit(1)
	leftfile = sys.argv[1]
	rightdir = sys.argv[2]
	extension = sys.argv[3]

	pattern = "%s/*.%s" % (rightdir, extension)
	
	leftname = os.path.basename(leftfile)

	rightfiles = glob(pattern)
	rightmap = {}
	for rightfile in rightfiles:
		rightname = os.path.basename(rightfile)
		if rightname == leftname:
			continue
		rightmap[rightname] = rightfile

# 	leftfile = "testinput/csp6i.mm9.bed"
# 	rightfile = "testinput/random10.bed"
# #	rightmap = {	"random10"	:	"testinput/random10.bed"	}
# 	rightmap = {	"H2B"	:	"testinput/H2B-pileup-wig.txt"	}
# 	
 	print "working on right files", rightmap.keys()

	# force it not to sort right, to test bad sort reporting
	sortedfiles = [ "H2B" ]
	og = overlapBedFilesPlaneSweepMany(leftfile, rightmap, sortedfiles=sortedfiles)
#	og = overlapBedFilesPlaneSweep(leftfile, rightfile)

	outdir = "testoutput"
	fhs = {}
	for name in rightmap:
		outfile = "%s/ogtester-%s.txt" % (outdir, name)
		ofh = open(outfile, "w")
		fhs[name] = ofh
	for overlap in og:
		rightname = overlap[1][0]
		ofh = fhs[rightname]
		print >> ofh, overlap


# 	outfile = "testoutput/ogtester-random10.bed.txt"
# 	ofh = open(outfile, "w")
# 	for overlap in og:
# 		print >> ofh, overlap
