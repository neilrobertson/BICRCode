#!/usr/bin/env python

import os

linkpath = "linked.py"
filepath = "combitester.py"

if os.access(linkpath, os.R_OK):
	os.remove(linkpath)
os.symlink(filepath, linkpath)
