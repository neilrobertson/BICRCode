#!/usr/bin/env python

import numpy as np
import csv
import resource

import array

def readAffy(fname, points, chrom):
	fh = open(fname, "r")
	reader = csv.reader(fh, delimiter="\t")
	reference_dict = {}
	count = 0
	for row in reader:
		chrom, centre, data = row
		centre = int(centre)
 		points.append(float(data))
		chrom = reference_dict.setdefault(chrom,chrom)
		coords[count] = (chrom,centre)
		count += 1
	fh.close()
	

def g2arr(fname):
	# the 'S100' should be modified to be large enough for your string field.
	dt = np.dtype({'names': ['chromo', 'position', 'dpoint'], 'formats': ['S6', 'i4', 'f8']})
	return np.loadtxt(fname, delimiter='\t', dtype=dt) 
	
if __name__ == '__main__':
	fname = "/home/pzs/histone/HISTONE_DATA/raw_affy_data/K562_H3K36me1_Normtogether_1_signal.txt-processed.txt"
#	fname = "tester.txt"
	# affy method
# 	arr = g2arr(fname)
# 	print len(arr['chromo'])
# 	print len(arr['position'])
# 	print len(arr['dpoint'])

	points = array.array('f')
	coords = {}
	readAffy(fname, points, coords)
	print resource.getrusage(resource.RUSAGE_SELF)
	raw_input()
