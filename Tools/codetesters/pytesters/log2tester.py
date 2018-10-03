#!/usr/bin/env python

from numpy import array
from scipy.stats import zs

import csv

reader = csv.reader(open("K562_H3K27me_mono-tri.txt", "rb"), delimiter="\t")
writer = csv.writer(open("K562_H3K27me_mono-tri.zscored.tester.txt", "wb"), delimiter="\t")

ids = []
vals = []
for row in reader:
	thisid, val = row
	ids.append(thisid)
	vals.append(float(val))
	
newvals = zs(array(vals))

for i in range(len(ids)):
	thisid = ids[i]
	thisval = newvals[i]
	writer.writerow((thisid, thisval))
	
	
