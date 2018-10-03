#!/usr/bin/env python

import csv

coordreader = csv.reader(open("889genes.csv", "r"), delimiter=" ")

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

genecoords = {}
for row in coordreader:
	if not row:
		continue
	chrom, start, end, gname = row[:4]
	gcoord = "%s:%s-%s" % (chrom, start, end)
	genecoords[gname] = gcoord
	
exonreader = csv.reader(open("/home/pzs/histone/HISTONE_DATA/data/microarray/Monocyte_EXONS.txt", "r"), delimiter="\t")
exonreader.next()

writer = csv.writer(open("humanexonsMonocyte.csv", "wb"), delimiter="\t")

for row in exonreader:
	if not row:
		continue
	try:
		gname = row[6]
		tname = row[6]
		ename = row[0]
		echrom = row[2]
		estart = row[3]
		eend = row[4]
	except:
		print row
		raise
	ecoord = "%s:%s-%s" % (echrom, estart, eend)
	tcoord = genecoords[gname]
	outrow = [ gname, tname, tcoord, ename, ecoord ]
	writer.writerow(outrow)

