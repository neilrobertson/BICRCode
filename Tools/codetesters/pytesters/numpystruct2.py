#!/usr/bin/env python

import numpy as np

def g2arr(fname):
	# the 'S100' should be modified to be large enough for your string field.
#	dt = np.dtype({'names': ['chromo', 'position', 'dpoint'], 'formats': ['S6', np.int, np.float]})
	dt = {'names': ['chromo', 'position', 'dpoint'], 'formats': ['S6', np.int, np.float]}
	fh = open(fname, "r")
	return np.loadtxt(fh, delimiter='\t', dtype=dt) 
	
if __name__ == '__main__':
#	fname = "/home/pzs/histone/HISTONE_DATA/raw_affy_data/K562_H3K36me1_Normtogether_1_signal.txt-processed.txt"
	fname = "/home/pzs/histone/HISTONE_DATA/raw_affy_data/tester.txt"
	arr = g2arr(fname)
	print arr
	print arr['chromo']
	print arr['position']
	print arr['dpoint'] 
