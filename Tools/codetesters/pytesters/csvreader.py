#!/usr/bin/env python

import csv

filename = "/home/pzs/histone/selectedsamples/2210_02_K4Me1.csv";


reader = csv.reader(open(filename, "rb"))
for row in reader:
	try:
		print row[4]
	except:
		print "uh?"
