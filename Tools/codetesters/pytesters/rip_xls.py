#!/usr/bin/env python

import xlrd
import sys
import os

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

def dumpOneTab(tab, name):
	fh = open(name, "w")
	for rownum in range(tab.nrows):
		line = []
		for colnum in range(tab.ncols):
			cellval = tab.cell_value(rowx=rownum, colx=colnum)
			try:
				cellval = str(cellval)
			except:
				cellval = "--garbled--"
			line += [cellval]
		print >> fh, ",".join(line)
	fh.close()
			

if len(sys.argv) != 2:
	print "usage: ./rip_xls.py <file.xls>"
	sys.exit(1)

filename = sys.argv[1]
filebase = os.path.basename(filename)[:-4]

print "reading workbook..."
xls = xlrd.open_workbook(filename)

for sheetname in xls.sheet_names():
	dumpfile = "%s-%s.csv" % (filebase, sheetname)
	print "dumping tab", dumpfile, "..."
	sheet = xls.sheet_by_name(sheetname)
	dumpOneTab(sheet, dumpfile)

