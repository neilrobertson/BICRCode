#!/usr/bin/env python

from sets import Set
import csv

from breakpoint import *

leftfile = "/home/pzs/histone/HISTONE_DATA/data/microarray/K562_GENES.txt"
rightfile = "/home/pzs/histone/IAN_newbuild.csv"

leftindex = 0
rightindex = 0 

leftreader = csv.reader(open(leftfile, "rb"), delimiter="\t")
rightreader = csv.reader(open(rightfile, "rb"), delimiter=",")

leftset = Set()
rightset = Set()

for row in leftreader:
	leftset.add(row[leftindex])

for row in rightreader:
	rightset.add(row[rightindex])

breakpoint()

# for name in rightset:
# 	if name not in leftset:
# 		print name
