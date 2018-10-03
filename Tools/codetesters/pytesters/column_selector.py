#!/usr/bin/env python

# for extracing specific columns from csv files. Very simple, but we still need it.

# no doubt this will subject to a number of custom hacks for specific tasks

# try to keep the basic file column_selector.py intact and use custom_selector.py for hacks

import sys
import os.path
import getopt
from sets import Set

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

def processCSV(filename, outfile, columns, delimiter, append="", header=False, unique=None):
	import csv
	reader = csv.reader(open(filename, "rb"), delimiter=delimiter)
	if outfile:
		outfh = open(outfile, "wb")
	else:
		outfh = sys.stdout
	writer = csv.writer(outfh)
	if header:
		reader.next()
	if unique != None:
		seen = Set()
	for row in reader:
		if not row:
			continue
		if unique != None:
			if row[unique] in seen:
				continue
		try:
			outpoints = [ row[i] for i in columns ]
		except IndexError:
			print "could not index row:"
			print row
		if append:
			outpoints.append(append)
		writer.writerow(outpoints)
		if unique != None:
			seen.add(row[unique])
	
def processExcel(filename, outfile, columns):
	import xlrd
	book = xlrd.open_workbook(filename)
	writer = csv.writer(open(outfile, "wb"))
	if book.nsheets > 1:
		print "can't process excel file with more than one sheet"
		sys.exit(1)
	sheet = book.sheet_by_index(0)
	for rownum in xrange(sheet.nrows):
		writer.writerow([ sheet.cell_value(rownum, col) for col in columns ])
		

def usage():
	print "usage: ./column_selector.py [OPTIONS] <filename>"
	print "\t-h\t\t\tPrint this message"
	print "\t-t, --type=TYPE\t\ttype of input file {csv, tsv, ssv, xls}"
	print "\t-a, --append=STRING\titem to add at the end of each row"
	print "\t-o, --outfile=PATH\toutput file"
	print "\t-c, --columns=C1,C2,...\tcomma separated list of output columns"
	print "\t-d, --header\t\tstrip first line from input file"
	print "\t-u, --unique\t\tcolumn that we force to only contain one of each entry"
	sys.exit(1)
	

if __name__ == "__main__":

	if len(sys.argv) < 3:
		usage()
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ht:a:f:c:du:", ["help", "type=", "append=", "infile=", "columns=", "unique=", "header"])
	except getopt.GetoptError, err:
		# print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage()
		sys.exit(2)
	outfile = None
	append = ""
	unique = None
	delimiter = None
	header = False
	intype = None
	for o, a in opts:
		if o in ("-t", "--type"):
			intype = a
		elif o in ("-a", "--append"):
			append = a
		elif o in ("-o", "--outfile"):
			outfile = a
		elif o in ("-c", "--columns"):
			columns = a.split(",")
		elif o in ("-d", "--header"):
			header = True
		elif o in ("-u", "--unique"):
			try:
				unique = int(a)
			except ValueError:
				print "unique column must be an integer"
				sys.exit(1)
		elif o in ("-h", "--help"):
			usage()
			sys.exit()
		else:
			assert False, "unhandled option"

	if len(args) > 1:
		print "too many arguments"
		usage()

	filename = args[0]

	base, ext = os.path.splitext(filename)
	if not intype:
		# if no type specified, try to guess
		intype = ext[1:]

	if not os.access(filename, os.R_OK):
		print "could not find filename:", filename
		sys.exit(1)

# 	if not outfile:
# 		outfile = "%s%s.out.csv" % (base, "".join(columns))
	if outfile:
		print "dumping into file", outfile

	columns = map(int, columns)

	if intype.endswith("sv"):
		if intype[0] == "c":
			delimiter = ","
		elif intype[0] == "t":
			delimiter = "\t"
		elif intype[0] == "s":
			delimiter = " "
		processCSV(filename, outfile, columns, delimiter=delimiter, append=append, header=header, unique=unique)
	elif intype == "xls":
		processExcel(filename, outfile, columns)
	else:
		print "can't deal with file of type:", intype
		sys.exit(1)
