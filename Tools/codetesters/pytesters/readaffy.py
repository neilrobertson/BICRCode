#!/usr/bin/env python

import numpy
import csv
import cStringIO
import glob

from breakpoint import *

chrommap = {	'1'	:	"chr1", 
				'2'	:	"chr2", 
				'3'	:	"chr3", 
				'4'	:	"chr4", 
				'5'	:	"chr5", 
				'6'	:	"chr6", 
				'7'	:	"chr7", 
				'8'	:	"chr8", 
				'9'	:	"chr9", 
				'0'	:	"chr10", 
				'a'	:	"chr11", 
				'b'	:	"chr12", 
				'c'	:	"chr13", 
				'd'	:	"chr14", 
				'e'	:	"chr15", 
				'f'	:	"chr16", 
				'g'	:	"chr17", 
				'h'	:	"chr18", 
				'i'	:	"chr19", 
				'j'	:	"chr20", 
				'k'	:	"chr21", 
				'l'	:	"chr22", 
				'x'	:	"chrX", 
				'y'	:	"chrY"
	}				



def reverseMap(cmap):
	rmap = {}
	for key in chrommap:
		value = chrommap[key]
		rmap[value] = key
	return rmap

def getNumPoints(fh):
	wholefile = fh.read()
	fh.seek(0)
	return wholefile.count("\n")
	
def readAffyFile(fname):
	print "reading affy file", fname
	fh = open(fname)
	n = getNumPoints(fh)

	chromio = cStringIO.StringIO()
	coords = numpy.zeros(n, dtype=int)
	points = numpy.zeros(n)

	reader = csv.reader(fh, delimiter="\t")
	count = 0
	for row in reader:
		if not row:
			continue
		chrom, coord, point = row
		chromio.write(rmap[chrom])
		coords[count] = coord
		points[count] = point
		count += 1
	return (chromio, coords, points)

def readAffyFiles(fmap):
	pointsets = {}
	for treatment in fmap:
		fglob = fmap[treatment]
		chromios = []
		coordas = []
		pointas = []
		files = glob.glob(fglob)
		for fname in files:
			chromio, coords, points = readAffyFile(fname)
			chromios.append(chromio)
			coordas.append(coords)
			pointas.append(points)
		chroms = "".join([ cio.getvalue() for cio in chromios ])
		coords = numpy.concatenate([ x for x in coordas ])
		points = numpy.concatenate([ x for x in pointas ])
		pointsets[treatment] = (chroms, coords, points)
	return pointsets

rmap = reverseMap(chrommap)

basedir = "/home/pzs/histone/HISTONE_DATA/raw_affy_data/"
prefix = "K562_"
# includes glob *
suffix = "_Normtogether_*_signal.txt-processed.txt"
treatments = [ "H3K27me1", "H3K27me3", "H3K36me1", "H3K36me3" ]
filemap = {}
for treatment in treatments:
	filemap[treatment] = "%s%s%s%s" % (basedir, prefix, treatment, suffix)

psl = readAffyFiles(filemap)
breakpoint()
