#!/usr/bin/env python

import csv
from sets import Set
import sys

if len(sys.argv) != 3:
	print "usage: ./uniquecolumn.py <filename> <column>"
	sys.exit(1)
filename = sys.argv[1]
try:
	column = int(sys.argv[2])
except ValueError:
	print "column must be an integer"
	sys.exit(1)

reader = csv.reader(open(filename, "rb"), delimiter="\t")

ids = Set()
for row in reader:
	ids.add(row[column])
	
idslist = list(ids)
idslist.sort()
print "\n".join(idslist)
