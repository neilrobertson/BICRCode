#!/usr/bin/env python

import csv
from sets import Set
from breakpoint import *

expressionout = "K562exp.csv"
coordout = "K562coord.csv"

expcsv = csv.writer(open(expressionout, "w"))
coordout = csv.writer(open(coordout, "w"))


exonin = "/home/pzs/codetesters/pytesters/humanexonsK562.csv"
genein = "/home/pzs/histone/HISTONE_DATA/data/final/K562_EXPRESSION.csv"

exoncsv = csv.reader(open(exonin, "r"), delimiter="\t")
genecsv = csv.reader(open(genein, "r"))

# keep the geneids that we have exon coordinates for
geneids = Set()

for row in exoncsv:
	geneid = row[0]
	geneids.add(geneid)

genecsv.next()
	
for row in genecsv:
	geneid, alias, region, chrom, start, end, strand, activation = row
	if geneid in geneids and alias in geneids:
		print "warning: both geneid and alias in exon list", geneid, alias
	if geneid not in geneids:
		if alias in geneids:
			geneid = alias
		else:
			print "no exon coords for:", geneid, alias
			breakpoint()
	exprow = [ geneid, activation ]
	expcsv.writerow(exprow)
	end = int(end) + 1
	coordrow = [ geneid, chrom, start, end, strand ]
	coordout.writerow(coordrow)
