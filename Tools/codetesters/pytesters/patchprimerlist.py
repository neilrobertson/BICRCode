#!/usr/bin/env python

import csv
from Bio import SeqIO

from breakpoint import *

primerfile = "/home/pzs/primerdesign/primerdesign/tags/parallel/promoterprimers.csv"
outfile = "processedprimers.csv"

reader = csv.reader(open(primerfile, "r"))
writer = csv.writer(open(outfile, "w"))

for row in reader:
	rowlen = len(row)
	if rowlen == 4:
		assert(row[-1] == "site not present!")
		writer.writerow(row)
	elif rowlen == 6:
		assert(row[-1] == "None found!")
		writer.writerow(row)
	elif rowlen == 11:
		writer.writerow(row)
		continue
	elif rowlen == 10:
		left = row[5]
		right = SeqIO.Seq(row[6])
		right = str(right.reverse_complement())
		fullseq = row[4]
		leftindex = fullseq.index(left)
		rightindex = fullseq.index(right) + len(right)
		product = fullseq[leftindex:rightindex]
		row.insert(7, product)
		writer.writerow(row)
	else:
		print "unknown row type", row
