#!/usr/bin/env python

import csv
from breakpoint import *

exonfile = "/home/pzs/histone/HISTONE_DATA/uniq_exons_and_introns_NCBI35.txt"

nameindex = 0
idindex = None
strandindex = None
chromindex = 3
startindex = 4
endindex = 5


ereader = csv.reader(open(exonfile, "r"), delimiter="\t")

ereader.next()


ecoords = {}
gdetails = {}
lastname = ""
# read in all the exons
for row in ereader:
	gname = row[nameindex]
	chrom = row[chromindex]
	start = row[startindex]
	end = row[endindex]
	if gname != lastname and gname in ecoords:
		pass
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



exonfile = "/home/pzs/histone/HISTONE_DATA/EXONS_with_ENS_ID.txt"

nameindex = 6
idindex = 10
strandindex = 5
chromindex = 2
startindex = 3
endindex = 4


ereader = csv.reader(open(exonfile, "r"), delimiter="\t")

ereader.next()

ecoords2 = {}
gdetails2 = {}
lastname = ""
# read in all the exons
for row in ereader:
	gname = row[nameindex]
	chrom = row[chromindex]
	start = row[startindex]
	end = row[endindex]
	if gname != lastname and gname in ecoords2:
		pass
	if idindex != None:
		gid = row[idindex]
	if strandindex != None:
		strand = row[strandindex]
	thisecoords2 = ecoords2.setdefault(gname, [])
	thisecoords2.append((chrom, int(start), int(end)))
	details = [gname]
	if idindex != None:
		details.append(gid)
	if strandindex != None:
		details.append(strand)
	gdetails2[gname] = details
	lastname = gname


for gname in ecoords2:
	numexons2 = len(ecoords2[gname])
	gid = gdetails2[gname][1]
	if gid in ecoords:
		numexons = len(ecoords[gid])
		if numexons != numexons2:
			print "Exon discrepancy:", gid, "In Rob's list:", numexons2, "In 10000 list", numexons
