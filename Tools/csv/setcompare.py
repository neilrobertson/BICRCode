#!/usr/bin/env python

# shell utility to compare the columns from two csv files

import sys
import csv

from bp import breakpoint

def getOneSet(filename, column):
	fh = open(filename)
	try:
		# if the read size is too small, it can't identify the delimiter properly
		dialect = csv.Sniffer().sniff(fh.read(5000), delimiters=" ,\t")
	except csv.Error:
		# couldn't find a dialect. Maybe this is only one column (not technically csv)
		fh.seek(0)
		assert column == 0
		return set([ line.strip() for line in fh ])
	fh.seek(0)
	reader = csv.reader(open(filename), dialect=dialect)
	return set([ row[column] for row in reader if row ])

def summarise(set1, set1name, set2, set2name):
	print "%s has %d items" % (set1name, len(set1))
	print "%s has %d items" % (set2name, len(set2))
	print "present in %s but not in %s: %d" % (set1name, set2name, len(set1 - set2))
	print "present in %s but not in %s: %d" % (set2name, set1name, len(set2 - set1))
	print "intersection has %d items" % (len(set1.intersection(set2)))
	print "union has %d items" % (len(set1.union(set2)))

if __name__ == "__main__":
	outputtypes = [ "summary", "intersection", "union", "leftonly", "rightonly" ]
	if len(sys.argv) != 6 or sys.argv[5] not in outputtypes:
		print "usage: setcompare.py <file1> <f1-column> <file2> <f2-column> <outputtype>"
		print "output types:", outputtypes
		sys.exit(1)
	set1name = sys.argv[1]
	set1col = int(sys.argv[2])
	set2name = sys.argv[3]
	set2col = int(sys.argv[4])
	set1 = getOneSet(set1name, set1col)
	set2 = getOneSet(set2name, set2col)
	outputtype = sys.argv[5]
	if outputtype == "summary":
		summarise(set1, set1name, set2, set2name)
	elif outputtype == "intersection":
		inter = set1.intersection(set2)
		for item in inter:
			print item
	elif outputtype == "union":
		union = set1.union(set2)
		for item in union:
			print item
	elif outputtype == "leftonly":
		leftonly = set1 - set2
		for item in leftonly:
			print item
	elif outputtype == "rightonly":
		rightonly = set2 - set1
		for item in rightonly:
			print item
	else:
		print "valid analysis types are:", analysistypes
