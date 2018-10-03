#!/usr/bin/env python

import csv
from breakpoint import *

exonfile = "newexonlist.csv"
genefile = "/home/pzs/histone/HISTONE_DATA/IAN_GENES.csv"

outfile = "IAN_GENES_butchered.csv"

ereader = csv.reader(open(exonfile, "r"), delimiter="\t")
greader = csv.reader(open(genefile, "r"))

bwriter = csv.writer(open(outfile, "w"))

ecoords = {}
gcoords = {}
# read in all the exons
for row in ereader:
	gid = row[0]
	ecoord = row[4]
	chrom, numbers = ecoord.split(":")
	start, end = numbers.split("-")
	thisecoords = ecoords.setdefault(gid, [])
	thisecoords.append((chrom, int(start), int(end)))


# merge the exons together to form a complete extent
for genename in ecoords:
	thisecoord = ecoords[genename]
	# take a sample as a start
	finalchrom, finalstart, finalend = thisecoord[0]
	for ecoord in thisecoord:
		chrom, start, end = ecoord
		assert(chrom == finalchrom)
		finalstart = min(finalstart, start)
		finalend = max(finalend, end)
	gcoords[genename] = (finalchrom, finalstart, finalend)

# for each gene listed
for row in greader:	
	name, chrom, start, end, strand, geneid = row
	start = int(start)
	end = int(end)
	if geneid in gcoords:
		echrom, estart, eend = gcoords[geneid]
	else:
		if name in gcoords:
			echrom, estart, eend = gcoords[name]
		else:
			print "no entry for", geneid, name
	# if the extent is bad, change it
	try:
		assert(echrom == chrom)
	except AssertionError:
		print "non-matching chromosome for gene", name, geneid
	if start != estart or end != eend:
		print "amending gene", geneid, name
		newrow = [ name, echrom, estart, eend, strand, geneid ]
		bwriter.writerow(newrow)
	else:
		bwriter.writerow(row)
