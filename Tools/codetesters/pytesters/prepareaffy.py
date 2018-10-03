#!/usr/bin/env python

import csv

basedir = "/home/pzs/histone/HISTONE_DATA/raw_affy_data/"
infiles = [ "K562_H3K27me1_Normtogether_1_signal.txt", "K562_H3K27me3_Normtogether_4_signal.txt", "K562_H3K36me1_Normtogether_7_signal.txt", "K562_H3K27me1_Normtogether_2_signal.txt", "K562_H3K27me3_Normtogether_5_signal.txt", "K562_H3K36me3_Normtogether_1_signal.txt", "K562_H3K27me1_Normtogether_3_signal.txt", "K562_H3K27me3_Normtogether_6_signal.txt", "K562_H3K36me3_Normtogether_2_signal.txt", "K562_H3K27me1_Normtogether_4_signal.txt", "K562_H3K27me3_Normtogether_7_signal.txt", "K562_H3K36me3_Normtogether_3_signal.txt", "K562_H3K27me1_Normtogether_5_signal.txt", "K562_H3K36me1_Normtogether_1_signal.txt", "K562_H3K36me3_Normtogether_4_signal.txt", "K562_H3K27me1_Normtogether_6_signal.txt", "K562_H3K36me1_Normtogether_2_signal.txt", "K562_H3K36me3_Normtogether_5_signal.txt", "K562_H3K27me1_Normtogether_7_signal.txt", "K562_H3K36me1_Normtogether_3_signal.txt", "K562_H3K36me3_Normtogether_6_signal.txt", "K562_H3K27me3_Normtogether_1_signal.txt", "K562_H3K36me1_Normtogether_4_signal.txt", "K562_H3K36me3_Normtogether_7_signal.txt", "K562_H3K27me3_Normtogether_2_signal.txt", "K562_H3K36me1_Normtogether_5_signal.txt", "K562_H3K27me3_Normtogether_3_signal.txt", "K562_H3K36me1_Normtogether_6_signal.txt" ]


for afile in infiles:
	apath = basedir + afile
	print "working at file", apath
	abase = apath[:-4]
	outpath = apath + "-processed.txt"
	ifh = open(apath, "r")
	ofh = open(outpath, "w")
	for line in ifh:
		if line.startswith("chr"):
			print >> ofh, line[:-1]
	ifh.close()
	ofh.close()
	
		
