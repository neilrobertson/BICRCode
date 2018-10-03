#!/usr/bin/env python

import csv

infile = "/home/pzs/codetesters/perltesters/gccounts.csv"

reader = csv.reader(open(infile, "r"))

for row in reader:
	try:
		this_id, seq, count = row
		gcount = seq.count("G")
		ccount = seq.count("C")
		assert(gcount + ccount == int(count))
	except:
		print row
		print gcount
		print ccount
		print count
		raise
