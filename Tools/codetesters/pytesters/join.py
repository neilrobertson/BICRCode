#!/usr/bin/env python

# another utility for doing SQL like things on csv files.
# in this case, look up a list of IDs in one file
# dump the lines from the other file based on these ids

import csv
from sets import Set

fromfile = "/home/pzs/histone/HISTONE_DATA/data/annotation/IAN_GENES.txt"
tofile = "/home/pzs/histone/HISTONE_DATA/K562_GENES.25.txt"
outfile = "K562_genes_618.txt"

joinfrom = Set()

reader1 = csv.reader(open(fromfile, "rb"), delimiter=" ")
for row in reader1:
	if not row:
		continue
	thisid = row[3]
	joinfrom.add(thisid)
	
writer = csv.writer(open(outfile, "wb"), delimiter="\t")

reader2 = csv.reader(open(tofile, "rb"), delimiter="\t")
for row in reader2:
	if row[0] in joinfrom:
		writer.writerow(row)
		

