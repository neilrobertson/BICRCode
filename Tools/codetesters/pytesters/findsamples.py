#!/usr/bin/env python

import csv
import os
from shutil import copy

def stripName(name):
	removes = ["OK", "Weak"]
	for remove in removes:
		name = name.replace(remove, "")
	name = name.replace(" ", "")
	name = name.strip("_")
	return name

basedirectory = "/home/pzs/histone/alldata"
directories = [ "BIOREP1", "BIOREP2", "BIOREP3", "BIOREP4", "BIOREP5" ] 
targetdirectory = "/home/pzs/histone/selectedsamples"
samplelist = targetdirectory + "/selectedsamples.csv"

# parse to find the samples we want
reader = csv.reader(open(samplelist, "rb"))
filenames = []
for row in reader:
	if not row[1]:	
		continue
	antibody = row[0]
	for exp in row[1:]:
		filename = stripName(exp) + "_" + antibody + ".csv"
		filenames.append(filename)

fileindex = {}
# index file space
for directory in directories:
	thispath = basedirectory + "/" + directory
	listing = os.listdir(thispath)
	for filename in listing:
		filepath = thispath + "/" + filename
		fileindex[filename] = filepath

# look for those samples in the different directory
for filename in filenames:
	if filename in fileindex:
		frompath = fileindex[filename]
		topath = targetdirectory + "/" + filename
		print "found file:", filename, "copying to path", topath
		copy(frompath, topath)
	else:
		print "could not find file:", filename
	
