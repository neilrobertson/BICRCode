#!/usr/bin/env python

import csv
import sys

exonlist = "/home/pzs/histone/HISTONE_DATA/EXONS_with_ENS_ID.txt"

reader = csv.reader(open(exonlist, "r"), delimiter="\t")
writer = csv.writer(sys.stdout, delimiter=",")

genelist = [ "PDZK1IP1" ]

allgenes = {}

for row in reader:
	chrom = row[2]
	start = row[3]
	end = row[4]
	strand = row[5]
	name = row[6]
	geneid = row[10]
	allgenes[name] = [ name, chrom, start, end, strand, geneid ]
	
for gene in genelist:
	writer.writerow(allgenes[gene])
