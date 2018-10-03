#!/usr/bin/env python

import glob
import csv

basedir = "/home/pzs/publicdata/zhao09HDAC/"
gpatt = "*vstep.txt"


files = glob.glob(basedir + gpatt)

for vfile in files:
	print "working at file", vfile
	outfile = vfile + "-processed.txt"
	reader = csv.reader(open(vfile), delimiter="\t")
	writer = csv.writer(open(outfile, "w"), delimiter="\t")
	chrom = None
	for row in reader:
		if len(row) == 1:
			line = row[0]
			if "chrom" in line:
				pieces = line.split(" ")
				chrom = pieces[1].split("=")[1]
			else:
				continue
		elif len(row) == 2:
			assert chrom != None
			writer.writerow([ chrom ] + row)
