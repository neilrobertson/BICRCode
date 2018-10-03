#!/usr/bin/env python

import csv
from sets import Set

outfile = "/home/pzs/primerdesign/primerdesignpython/trunk/promoterprimers.csv"

print "looking up ids we've already found"
seen = Set()
prepreader = csv.reader(open(outfile, "rb"))
# strip off header
prepreader.next()
for line in prepreader:
	thisid = line[0]
	seen.add(thisid)

print seen
