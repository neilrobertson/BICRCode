#!/usr/bin/env python

import os

inpath = "/home/pzs/histone/compositeprofile/trunk/compositeprofile"

def oneDir(accumulator, dirname, files):
	for file_ in files:
		if file_.endswith(".py"):
			accumulator.append("/".join([ dirname, file_ ]))
		
dirlist = []
os.path.walk(inpath, oneDir, dirlist)
print dirlist

