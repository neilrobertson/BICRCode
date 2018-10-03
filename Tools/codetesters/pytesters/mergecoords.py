#!/usr/bin/env python

import csv

filename = "/home/pzs/histone/HISTONE_DATA/data/annotation/ENCODE_union_annotation.txt"
outfile = "/home/pzs/histone/human_ids.txt"

reader = csv.reader(open(filename, "rb"), delimiter="\t")
writer = csv.writer(open(outfile, "wb"), delimiter="\t")

for row in reader:
	chrom = row[2]
	start = row[3]
	end = row[4]
	coord = "chr%s:%s-%s" % (chrom, start, end)
	newrow = [ row[0], coord ]
	writer.writerow(newrow)
