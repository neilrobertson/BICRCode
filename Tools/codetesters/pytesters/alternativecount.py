#!/usr/bin/env python

import sys
import csv
from pylab import *

from breakpoint import *

sys.path.append("/home/pzs/histone/compositeprofile/trunk")

from bins import AbstractSpreadBin
from misc import doesOverlap

infile = "/home/pzs/histone/HISTONE_DATA/uniq_exons_and_introns_NCBI35.txt"
numbins = 20


def buildBins(numbins):
	bins = []
	for i in range(numbins):
		b = AbstractSpreadBin(i, i+1)
		bins.append(b)
	leftbin = AbstractSpreadBin(-1, 0)
	rightbin = AbstractSpreadBin(numbins, numbins + 1)
	return bins

def fileIntoBins(bins, point):
	start, end, data, pointtype = point
	for bin in bins:
		if doesOverlap(start, end, bin.lower, bin.upper):
			bin.drop(data, pointtype)

def updatePositions(genepoints):
	# find highest numbered exon and intron
	maxi = -1
	maxe = -1
	for point in genepoints:
		position, data, pointtype = point
		if position != "last":
			if pointtype == "canonical" or pointtype == "alternative":
				maxe = max(maxe, int(position))
			elif pointtype == "intron":
				maxi = max(maxi, int(position))
	# add one for the "last" that we haven't yet counted
	maxe += 1
	# normalise them around the total number of exons or introns
	newpoints = []
	for point in genepoints:
		position, data, pointtype = point
		if position == "last":
			if pointtype == "canonical" or pointtype == "alternative":
				start = float(maxe - 1) / maxe
			elif pointtype == "intron":
				start = float(maxi - 1) / maxi
			end = 1
		else:
			if pointtype == "canonical" or pointtype == "alternative":
				end = float(position) / maxe
				start = (float(position) - 1.0) / maxe
			elif pointtype == "intron":
				end = float(position) / maxi
				start = (float(position) - 1.0) / maxi
		start *= numbins
		end *= numbins
		newpoints.append((start, end, data, pointtype))
	return newpoints

def binplot(bins):
	xs = []
	introns = []
	canonical = []
	alternative = []
	for bin in bins:
		xs.append(bin.midx)
		introns.append(len(bin.getRawData("intron")))
		canonical.append(len(bin.getRawData("canonical")))
		alternative.append(len(bin.getRawData("alternative")))

	figure(1)

	plot(xs, introns)
	plot(xs, canonical)
	plot(xs, alternative)
	legend(("introns", "canonical", "alternative"))

	show()

bins = buildBins(numbins)

reader = csv.reader(open(infile, "r"), delimiter="\t")

genepoints = []
lastgenename = ""
for row in reader:
	genename = row[0]
	# strip off E or I to give just the position
	position = row[1][1:]
	# exon. Strip off first character to get canonical/alternative designation
	if row[2].startswith("E"):
		eitype = row[2][1:]
	else:
		eitype = "intron"
	# finished reading in a whole gene. File points
	if lastgenename and lastgenename != genename:
		genepoints = updatePositions(genepoints)
		for point in genepoints:
			fileIntoBins(bins, point)
		genepoints = []
		maxposition = -1
	genepoints.append((position, 1, eitype))
	lastgenename = genename
# last one
# genepoints = updatePositions(genepoints)
# for point in genepoints:
# 	fileIntoBins(bins, point)

binplot(bins)
