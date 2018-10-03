#!/usr/bin/env python

import csv

infile = "/home/pzs/histone/wigglefiles/2210_02_K4Me1.wig.txt"

fh = open(infile, "rb")

line = fh.readline()
while not line.startswith("chr"):
	line = fh.readline()
	
dialect = csv.Sniffer().sniff(fh.read(1024))
fh.seek(0)
print dir(dialect)
reader = csv.reader(fh, dialect)

for line in reader:
	print line
