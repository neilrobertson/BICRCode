#!/usr/bin/env python

import csv

coordreader = csv.reader(open("889genes.csv", "rb"), delimiter=" ")

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
	
exonreader = csv.reader(open("interimexons.csv", "rb"))
exonreader.next()

writer = csv.writer(open("humanexonsENCODE.csv", "wb"), delimiter="\t")

for row in exonreader:
	if not row:
		continue
	gname, tname, ename, echrom, estart, eend = row
	ecoord = "%s:%s-%s" % (echrom, estart, eend)
	tcoord = genecoords[gname]
	outrow = [ gname, tname, tcoord, ename, ecoord ]
	writer.writerow(outrow)

