#!/usr/bin/env python

import os
import sys

host = os.uname()[1]
approot = os.path.abspath(os.path.dirname(__file__))
topdir = approot.split("/")[-1]
if host == "pzs-desktop": 
	if topdir == "trunk" or topdir == "pytesters":
		sys.path.append("/home/pzs/histone/compositeprofile/trunk")
	elif topdir == "deployed":
		sys.path.append("/home/pzs/histone/compositeprofile/deployed")
	weburl = "http://130.209.136.101/"
elif host == "mimas.nesc.gla.ac.uk":
	sys.path.append("/home/pzs/lib")
	weburl = "http://mimas.nesc.gla.ac.uk/~pzs/microarray/"

#### 
# settings
#####

apacheroot = "/var/www"
heroic_file = approot + "/datafiles/HEROIC_all.txt"
idfilesettings = { "coord" : 10, "id" : 5, "header" : True }
outputdir = "/files"
# use this one when we're outputting into the directory
localoutputdir = approot + outputdir
appwebroot = approot[len(apacheroot):]
# use this one when we're linking to it from a web page
weboutputdir = appwebroot + outputdir
