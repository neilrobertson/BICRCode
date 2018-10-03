#!/usr/bin/env python

import csv

filename = "/home/pzs/codetesters/perltesters/humanexonlist.csv"
outfilename = "/home/pzs/codetesters/pytesters/humanexonlist.csv"

reader = csv.reader(open(filename, "rb"), delimiter="\t")
writer = csv.writer(open(outfilename, "wb"), delimiter="\t")

for row in reader:
	outrow = row[:]
	outrow[2] = "chr" + row[2]
	outrow[4] = "chr" + row[4]
	writer.writerow(outrow)
