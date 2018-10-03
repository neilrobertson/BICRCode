#!/usr/bin/env python

import csv
from sets import Set

coordfile = "/home/pzs/codetesters/pytesters/promoterseqs.csv"
infile = "/home/pzs/primerdesign/primerdesign/trunk/promoterprimers.csv"
outfile = "promoterprimerslabelled.csv"

coordreader = csv.reader(open(coordfile, "r"))


coordmap = {}
for row in coordreader:
	coordmap[row[0]] = row[1]

reader = csv.reader(open(infile, "r"))

writer = csv.writer(open(outfile, "w"))

seen = Set()
reader.next()
for row in reader:
	second = row[1]
	if not second.startswith("chr"):
		thisid = row[0]
		if thisid in seen:
			continue
		seen.add(thisid)
		thiscoord = coordmap[thisid]
		row.insert(1, thiscoord)
	writer.writerow(row)

