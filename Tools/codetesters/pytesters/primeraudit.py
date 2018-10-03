#!/usr/bin/env python

import csv

# primerfile = "/home/pzs/primerdesign/primerdesign/tags/parallel/promoterprimers.csv"
primerfile = "processedprimers.csv"


reader = csv.reader(open(primerfile, "r"))

nositecount = 0
noprimerscount = 0
primerscount = 0
badrowcount = 0
noproductcount = 0
for row in reader:
	rowlen = len(row)
	if rowlen == 4:
		assert(row[-1] == "site not present!")
		nositecount += 1
	elif rowlen == 6:
		assert(row[-1] == "None found!")
		noprimerscount += 1
	elif rowlen == 11:
		primerscount += 1
	elif rowlen == 10:
		noproductcount += 1
	else:
		badrowcount += 1
		
print "no sites:", nositecount
print "no primers:", noprimerscount
print "primers:", primerscount
print "bad rows:", badrowcount
print "total:", nositecount + noprimerscount + primerscount + badrowcount
