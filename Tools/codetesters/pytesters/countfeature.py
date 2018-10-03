#!/usr/bin/env python

import csv

infile = "promoterprimerslabelled.csv"
outfile = "promoterprimers.csv"

reader = csv.reader(open(infile, "r"))
writer = csv.writer(open(outfile, "w"))

for row in reader:
	rowlen = len(row)
	if rowlen == 8:
		continue
	else:
		writer.writerow(row)
		
