#!/usr/bin/env python

import os

basedir = "/home/pzs/histone/selectedsamples"

allsamples = os.listdir(basedir)
allsamples.remove("selectedsamples.csv")

groupids = [ "H2B", "H3", "K4Me1", "K4Me3", "K9", "K27Me1", "K27Me3", "K36Me1", "K36Me3", "K79Me1", "K79Me3", "Goat", "R" ]
groups = {}
for filename in allsamples:
	ext = filename.split(".")[-1]
	if ext != "csv":
		continue
	for groupid in groupids:
		if groupid in filename:
			glist = groups.setdefault(groupid, [])
			glist.append(filename)

print groups
