#!/usr/bin/env python

import csv

infile = "/home/pzs/histone/HISTONE_DATA/data/annotation/IAN_GENES.txt"

chromindex = 0
startindex = 1
endindex = 2
nameindex = 3
strandindex = 5

reader = csv.reader(open(infile, "rb"), delimiter=" ")

for row in reader:
	if not row:
		continue
	try:
		chrom = row[chromindex]
		start = row[startindex]
		end = row[endindex]
		name = row[nameindex]
		strand = row[strandindex]
	except IndexError:
		print "bad line:", row
		continue
	try:
		int(start)
		int(end)
	except ValueError:
		# skip lines that don't have the right entries
		continue
	print "%s\t%s:%s-%s\t%s" % (name, chrom, start, end, strand)
