#!/usr/bin/env python

import os
import glob

basedir = "/home/pzs/histone/HISTONE_DATA/raw_affy_data/"
fpattern = basedir + "*-processed.txt"

allfiles = glob.glob(fpattern)


for fname in allfiles:
	newfile = fname.split(".")[0] + "-cropped.txt"
	fpath = fname
	print "working at path", fpath
	ifh = open(fpath, "r")
	ofh = open(newfile, "w")
	count = 0
	while count < 1000:
		line = ifh.readline()
		print >> ofh, line[:-1]
		count += 1
