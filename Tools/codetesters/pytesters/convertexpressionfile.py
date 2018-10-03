#!/usr/bin/env python

import csv
from sets import Set
from breakpoint import *

expressionout = "exp.csv"

expoutcsv = csv.writer(open(expressionout, "w"))


mapin = "/home/pzs/histone/HISTONE_DATA/IAN_GENES.csv"
expin = "/home/pzs/histone/HISTONE_DATA/data/final/Monocyte_EXPRESSION.csv"

mapcsv = csv.reader(open(mapin, "r"), delimiter=",")
expcsv = csv.reader(open(expin, "r"), delimiter=",")

# keep the geneids that we have exon coordinates for
geneidmap = {}

for row in mapcsv:
	genename = row[0]
	if not genename:
		continue
	try:
		geneid = row[5]
	except IndexError:
		geneid = ""
	if geneid:
		geneidmap[genename] = geneid
	else:
		geneidmap[genename] = genename

expcsv.next()
	
for row in expcsv:
	genename, alias, strand, activation = row
	if genename in geneidmap and alias in geneidmap:
		print "warning: both geneid and alias in geneidmap list", genename, alias
	if genename not in geneidmap:
		if alias in geneidmap:
			genename = alias
		else:
			print "no geneidmap coords for:", genename, alias
			continue
	geneid = geneidmap[genename]
	exprow = [ geneid, activation ]
	expoutcsv.writerow(exprow)
