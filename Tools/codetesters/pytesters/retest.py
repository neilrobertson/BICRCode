#!/usr/bin/env python

import re

infile = "/home/pzs/primerdesign/primer_design/1_amplicons/data/temp/primer3_output.temp"

fh = open(infile, "rb")

for line in fh:
	var, val = line.split("=")
	if re.match("PRIMER_RIGHT(_\d+)?$", var):
		print var
