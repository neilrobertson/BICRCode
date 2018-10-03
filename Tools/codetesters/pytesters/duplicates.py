#!/usr/bin/env python

import csv
from sets import Set

reader = csv.reader(open("/home/pzs/histone/HISTONE_DATA/IAN_GENES.csv", "r"))

allseen = Set()
for row in reader:
	thisid = row[5]
	if not thisid:
		continue
	if thisid in allseen:
		print "duplicate id", thisid
	else:
		allseen.add(thisid)
