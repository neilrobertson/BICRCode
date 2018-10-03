#!/usr/bin/env python

import csv

reader = csv.reader(open("/var/www/pythonapps/datafiles/HEROIC_all.txt", "rb"), delimiter="\t")

starts = []
reader.next()
for row in reader:
	coord = row[-1]
	try:
		chrom, pos = coord.split(":")
		start, end = pos.split("-")
	except ValueError:
		print "problem with line:", coord
	starts.append(int(start))
	
starts.sort()

diffs = []
for i in range(1, len(starts)):
	diff = starts[i] - starts[i-1]
	if diff < 10000:
		diffs.append(diff)
	
print sum(diffs)/len(diffs)
