#!/usr/bin/env python

# for extracing specific columns from csv files. Very simple, but we still need it.

# no doubt this will subject to a number of custom hacks for specific tasks

# try to keep the basic file column_selector.py intact and use custom_selector.py for hacks

import csv
import xlrd
import sys
import os.path

rfilter =  ['ENSMUSG00000052217', 'ENSMUSG00000052217', 'ENSMUSG00000052187', 'ENSMUSG00000050621', 'ENSMUSG00000019558', 'ENSMUSG00000073965', 'ENSMUSG00000045824', 'ENSMUSG00000073962', 'ENSMUSG00000073956', 'ENSMUSG00000073956', 'ENSMUSG00000073954', 'ENSMUSG00000061626', 'ENSMUSG00000051618', 'ENSMUSG00000073925', 'ENSMUSG00000060105', 'ENSMUSG00000016213', 'ENSMUSG00000052217', 'ENSMUSG00000055124', 'ENSMUSG00000073931', 'ENSMUSG00000056144', 'ENSMUSG00000056144', 'ENSMUSG00000079529']


def processCSV(filename, outfile, columns, delimiter):
	reader = csv.reader(open(filename, "rb"), delimiter=delimiter)
	writer = csv.writer(open(outfile, "wb"))
	for row in reader:
		outpoints = [ row[i] for i in columns ]
		writer.writerow(outpoints)
	
def processExcel(filename, outfile, columns):
	book = xlrd.open_workbook(filename)
	ofh = open(outfile, "w")
	if book.nsheets > 1:
		print "can't process excel file with more than one sheet"
		sys.exit(1)
	sheet = book.sheet_by_index(0)
	for rownum in xrange(sheet.nrows):
		row = [ sheet.cell_value(rownum, col) for col in columns ]
		if row[1] in rfilter:
			strand, ensid, chrom, start, end = row
			try:
				chrom = int(chrom)
			except ValueError:
				chrom = str(chrom)
			start = int(start)
			end = int(end)
			strand = int(strand)
			print >> ofh, "chr%s:%s-%s,%s\t%s" % (chrom, start, end, strand, ensid)
		

if __name__ == "__main__":

	if len(sys.argv) < 3:
		print "usage: ./column_selector.py <filename> <column_index1> <column_index2> ..."
		sys.exit(1)

	filename = sys.argv[1]
	if not os.access(filename, os.R_OK):
		print "could not find filename:", filename
		sys.exit(1)
	base, ext = os.path.splitext(filename)

	columns = sys.argv[2:]

	outfile = "missing.txt" 
	print "dumping into file", outfile

	columns = map(int, columns)

	if ext == ".csv":
		processCSV(filename, outfile, columns, ",")
	elif ext == ".tsv":
		processCSV(filename, outfile, columns, "\t")
	elif ext == ".xls":
		processExcel(filename, outfile, columns)
	else:
		print "can't deal with file of type:", ext
		sys.exit(1)
