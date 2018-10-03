#!/usr/bin/env python

import csv
from breakpoint import *

#exonfile = "/home/pzs/histone/HISTONE_DATA/EXONS_with_ENS_ID.txt"
exonfile = "/home/pzs/histone/HISTONE_DATA/uniq_exons_and_introns_NCBI35.txt"

nameindex = 0
idindex = None
strandindex = None
chromindex = 3
startindex = 4
endindex = 5


outfile = "forcedcoords10000.csv"

ereader = csv.reader(open(exonfile, "r"), delimiter="\t")

bwriter = csv.writer(open(outfile, "w"), delimiter="\t")

ereader.next()


ecoords = {}
gcoords = {}
gdetails = {}
lastname = ""
# read in all the exons
for row in ereader:
	gname = row[nameindex]
	chrom = row[chromindex]
	start = row[startindex]
	end = row[endindex]
	if gname != lastname and gname in ecoords:
		print "seen name before:", gname
		print chrom
	if idindex != None:
		gid = row[idindex]
	if strandindex != None:
		strand = row[strandindex]
	thisecoords = ecoords.setdefault(gname, [])
	thisecoords.append((chrom, int(start), int(end)))
	details = [gname]
	if idindex != None:
		details.append(gid)
	if strandindex != None:
		details.append(strand)
	gdetails[gname] = details
	lastname = gname


# merge the exons together to form a complete extent
for genename in ecoords:
	thisecoord = ecoords[genename]
	# take a sample as a start
	finalchrom, finalstart, finalend = thisecoord[0]
	for ecoord in thisecoord:
		chrom, start, end = ecoord
		try:
			assert(chrom == finalchrom)
		except AssertionError:
			finalchrom = "X/Y"
			finalstart = "???"
			finalend = "???"
			break
		finalstart = min(finalstart, start)
		finalend = max(finalend, end)
	gcoords[genename] = (finalchrom, finalstart, finalend)

for geneid in gcoords:
	row = gdetails[geneid] + list(gcoords[geneid])
	bwriter.writerow(row)
