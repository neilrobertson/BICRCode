#!/usr/bin/env python

import cPickle
import csv
import time
from sets import Set

from breakpoint import *
import intervaltree
import cintervaltree

usecintervaltree = False

def chromcmp(c1, c2):
	if c1==c2:
		return 0
	if c1 != c2:
		if c1.startswith("chr"):
			c1 = c1[3:]
		if c2.startswith("chr"):
			c2 = c2[3:]
		try:
			return cmp(int(c1), int(c2))
		except ValueError:
			return cmp(c1, c2)

def chromsort(x1, x2):
	l = x1[0]
	r = x2[0]
	if l != r:
		return chromcmp(l, r)
	else:
		return cmp(int(x1[1]), int(x2[1]))

def readTestSet(fname):
	points = []
	coords = []
	reader = csv.reader(open(fname, "r"), delimiter="\t")
	for row in reader:
		chrom, centre, point = row
		centre = int(centre)
		coords.append((chrom, centre-10, centre+10))
		points.append(point)
	return (coords, points)

def breakCoord(coord):
	try:
		(chrom,coords) = coord.split(":")
	except ValueError:
		print "failed to break coord:", coord
		breakpoint()
		assert(False)
	(start,end) = coords.split("-")
	return (chrom, int(start), int(end))

def initBins(genelist):
	maxbinsize = 0
	minbinsize = 99999999
	numbins = 0
	for genename in genelist:
		genecoord, bins = genelist[genename]
		newbins = []
		for bin in bins:
			binsize = bin[1] - bin[0]
			maxbinsize = max(maxbinsize, binsize)
			minbinsize = min(minbinsize, binsize)
			newbin = bin + ([],genename)
			newbins.append(newbin)
		numbins += len(bins)
		genelist[genename] = genecoord, newbins
	print "number of bins:", numbins
	print "max bin size:", maxbinsize
	print "min bin size:", minbinsize

def buildChromosomeList(genelist):
	# build reverse index to look up features by chromosome
	chromindex = {}
	for genename in genelist:
		genecoord, bins = genelist[genename]
		chrom, start, end = breakCoord(genecoord)
		chromlist = chromindex.setdefault(chrom, [])
		try:
			geneextent = "%s:%s-%s" % (chrom, bins[0][0], bins[-1][1])
		except IndexError:
			geneextent = genecoord
		chromlist.append((genename, geneextent, bins))
	for chrom in chromindex:
		chromlist = chromindex[chrom]
		chromlist.sort(chromsort)
	return chromindex

def buildIntervalTreesByChromosome(chromlist):
	print "building interval trees"
	its = {}
	t0 = time.clock()
	for chrom in chromlist:
		genelist = chromlist[chrom]
		allbins = []
		for gene in genelist:
			genename, genecoord, bins = gene
			allbins.extend(bins)
		if usecintervaltree:
			cit = cintervaltree.IntervalTree()
			for bin in allbins:
				cit.insert_interval(cintervaltree.Interval(bin[0], bin[1], value=(bin[2], bin[3])))
			its[chrom] = cit
		else:
			its[chrom] = intervaltree.IntervalTree(allbins)
	t1 = time.clock()
	print "tree contruction time: %.2f" % (t1 - t0)
	return its

def overlap(lb, ub, start, end):
	return not (ub < start or lb > end)

def firstBinIndex(lb, binstart, binsize):
	firstbinindex = int((lb - binstart) / binsize )
	if firstbinindex < 0:
		firstbinindex = 0
	return firstbinindex

def getFilingBinsByFirst(bins, lb, ub):
	fbins = []
	binstart = bins[0][0]
	binend = bins[0][1]
	firstbinindex = firstBinIndex(lb, binstart, binend)
	found = False
	for i in range(firstbinindex, len(bins)):
		bin = bins[i]
		binstart, binend, bincontents, genename = bin
		if overlap(lb, ub, binstart, binend):
			fbins.append(bin)
			found = True
		else:
			if found == True:
				break
	return fbins

def getFilingBinsNaive(bins, lb, ub):
	fbins = []
	found = False
	for bin in bins:
		binstart, binend, bincontents, genename = bin
		if overlap(lb, ub, binstart, binend):
			fbins.append(bin)
			found = True
		else:
			if found == True:
				break
	return fbins

def filePoint(flist, lb, ub, point, filecount):
	seen = False
	for feature in flist:
		featurename, fcoord, bins = feature
		chrom, start, end = breakCoord(fcoord)
		if overlap(lb, ub, start, end):
			seen = True
			bins = getFilingBinsByFirst(bins, lb, ub)
			for bin in bins:
				binstart, binend, bincontents, genename = bin
				bincontents.append(point)
			filecount.setdefault(featurename, 0)
			filecount[featurename] += len(bins)
		else:
			if seen:
				break

def filePointIntervalTree(ftree, lb, ub, point, filecount):
	seen = False
	bins = ftree.find(lb, ub)
	for bin in bins:
		if usecintervaltree:
			binstart, binend = bin.start, bin.end
			bincontents, genename = bin.value
		else:
			binstart, binend, bincontents, genename = bin
		bincontents.append(point)
		filecount.setdefault(genename, 0)
		filecount[genename] += len(bins)

def filePoints(coords, points, chromlist):
	filecount = {}
	for i in range(len(coords)):
		chrom, lb, ub = coords[i]
		point = points[i]
		flist = chromlist[chrom]
		filePoint(flist, lb, ub, point, filecount)
	for genename in filecount:
		print "genename", genename, "filed", filecount[genename], "times"

def filePointsIntervalTree(coords, points, itlist):
	filecount = {}
	for i in range(len(coords)):
		chrom, lb, ub = coords[i]
		point = points[i]
		ftree = itlist[chrom]
		filePointIntervalTree(ftree, lb, ub, point, filecount)
	for genename in filecount:
		print "genename", genename, "filed", filecount[genename], "times"

picfile = "/home/pzs/histone/compositeprofile/trunk/eistructure.pic"
testset = "tester2.txt"

fh = open(picfile, "r")

genelist = cPickle.load(fh)

initBins(genelist)

chromlist = buildChromosomeList(genelist)
its = buildIntervalTreesByChromosome(chromlist)

coords, points = readTestSet(testset)
t0 = time.clock()
#filePoints(coords, points, chromlist)
filePointsIntervalTree(coords, points, its)
t1 = time.clock()

print "filing time: %.2f" % (t1 - t0)
