#!/usr/bin/python

# parses a jacobian file, looking for the user-specified parameters 
# inserting parameter values from the prf supplied.
# alters the given file by writing into a temporary file and then copying
# does some error checking to make sure we see all the parameters in the prf


import sys
from xml.dom.minidom import parse
from sets import Set

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

parammarker = "#USERPARAMETERS"

# takes a prf file
# returns a dictionary mapping parameter name to value
def getParameters(prffile):
	params = {}
	dom = parse(prffile)
	paramelems = dom.getElementsByTagName("parameter")
	for paramelem in paramelems:
		name = paramelem.getAttribute("syntax_name")
		value = paramelem.getAttribute("value")
		params[name] = value
	return params
	

if len(sys.argv) != 3:
	print "usage: ./jacobian-prfinsert.py <prffile> <jacobianfile>"
	sys.exit(1)
		
prffile = sys.argv[1]
jacfile = sys.argv[2]
outfile = jacfile + ".out"

params = getParameters(prffile)
seen = Set()

insert = False

jfh = open(jacfile, "rb")
ofh = open(outfile, "wb")

count = 0
for line in jfh:
	count += 1
#	if count == 159:
#		breakpoint()
	# perl chop
	line = line.strip("\r\n")
	if not line or line.startswith("#") and line != parammarker:
		continue
	if insert:
		if line == parammarker:
			insert = False
			print >> ofh, line
			continue
		else:
			if "#" in line:
				line = line.split("#")[0]
			if ";" not in line:
				continue
			# get at the assignment
			try:
				assign, rest = line.split(";")
			except ValueError:
				print "problem splitting line"
				print line
				sys.exit(1)
			param, value = assign.split(":=")
			param = param.strip()
			value = value.strip()
			print "working at param", param
			if param not in params:
				print "parameter with no specified value:", param
			seen.add(param)
			try:
				newvalue = params[param]
			except KeyError:
				print "could not find parameter", param
				breakpoint()
			print >> ofh, "\t%s := %s;" % (param, params[param])
	if not insert:
		print >> ofh, line
		if line == parammarker:
			print "entering insert mode"
			if seen:
				print "More than one user parameter definition"
			insert = True

allparams = Set(params.keys())
if allparams != seen:
	print "Did not allocate all parameters"
	print allparams - seen
	
ofh.close()
