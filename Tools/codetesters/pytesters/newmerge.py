#!/usr/bin/env python

import csv
from breakpoint import *

outfile = "mergediscrepancies.csv"

exonfile = "/home/pzs/histone/HISTONE_DATA/uniq_exons_and_introns_NCBI35.txt"

nameindex = 0
idindex = None
strandindex = None
chromindex = 3
startindex = 4
endindex = 5
eiindex = 1

ereader = csv.reader(open(exonfile, "r"), delimiter="\t")

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
	ei = row[eiindex]
	if ei.startswith("I"):
		continue
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
		finalstart = min(finalstart, end)
		finalend = max(finalend, end)
		finalend = max(finalend, start)
	gcoords[genename] = (finalchrom, finalstart, finalend)

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
	thisecoords2 = ecoords2.setdefault(gid, [])
	thisecoords2.append((chrom, int(start), int(end)))
	lastname = gname


writer = csv.writer(open(outfile, "w"), delimiter="\t")

# from genes in Ian's list
genefile = "IAN_GENES_DAVE.csv"

greader = csv.reader(open(genefile, "r"))

# number of bases out the coordinates can be
tolerance = 100

writer.writerow(["genename", "geneid", "chrom", "start", "end", "classification", "coord match", "exon count match" ])
for row in greader:
	gname = row[0]
	gid = row[5]
	chrom = row[1]
	start = int(row[2])
	end = int(row[3])
	classification = row[8]
	outrow = [gname, gid, chrom, start, end, classification]
	# big file coordinate
	if gid in gcoords:
		bigchrom, bigstart, bigend = gcoords[gid]
		bigcoord = "%s:%s-%s" % (bigchrom, bigstart, bigend)
		coord = "%s:%s-%s" % (chrom, start, end)
		startdiff = abs(start - bigstart)
		enddiff = abs(end - bigend)
		if bigchrom != chrom:
			msg = "chromosomes do not match! 10k file: %s / IAN file: %s" % (bigcoord, coord)
			outrow.append(msg)
		elif startdiff > 100 or enddiff > 100:
			msg = "coordinates do not match! 10k file: %s / IAN file: %s" % (bigcoord, coord)
			outrow.append(msg)
		else:
			outrow.append("MATCH")
	else:
		outrow.append("geneid not in 10000 file")
	# number in Rob's file
	if gid in ecoords2:
		numexons2 = len(ecoords2[gid])
		if gid in ecoords:
			numexons = len(ecoords[gid])
			if numexons != numexons2:
				msg = "exon discrepancy: Rob's file: %s / 10k file: %s" % (numexons2, numexons)
				outrow.append(msg)
			else:
				outrow.append("MATCH")
		else:
			outrow.append("geneid not in 10000 file")
	else:
		outrow.append("geneid not in Rob's file")
	writer.writerow(outrow)
