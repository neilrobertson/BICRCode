#!/usr/bin/env python

import csv
from sets import Set

idcsv = "/home/pzs/histone/HISTONE_DATA/IAN_GENES_butchered.csv"
wholegenomeeilist = "/home/pzs/histone/HISTONE_DATA/uniq_exons_and_introns_NCBI35.txt"
oldeilist = "/home/pzs/histone/HISTONE_DATA/EXONS_with_ENS_ID.txt"

# We want a list where each row is:
# [ geneid, transcriptid, transcriptcoord, exonid, exoncoord, canonical/alternative]

# build a map from id to row that describes that row
# each row is [ id, id, genecoord, exonid, exoncoord, canonical/alternative ]
# the details from the whole genome list are preferred, so we read those first

geireader = csv.reader(open(wholegenomeeilist, "r"), delimiter="\t")
oldeireader = csv.reader(open(oldeilist, "r"), delimiter="\t")

ianreader = csv.reader(open(idcsv, "r"), delimiter=",")

genetorows = {}
genenametorows = {}
for row in geireader:
	ensid, ei, canalt, chrom, start, end = row
	# skip introns
	if ei.startswith("I"):
		continue
	coord = "%s:%s-%s" % (chrom, start, end)
	generows = genetorows.setdefault(ensid, [])
	generows.append([ ensid, ensid, ei, coord, canalt ])

for row in oldeireader:
	exonid = row[0]
	geneid = row[10]
	exonchrom = row[2]
	exonstart = row[3]
	exonend = row[4]
	genename = row[6]
	coord = "%s:%s-%s" % (exonchrom, exonstart, exonend)
	if geneid not in genetorows:
		generows = genetorows.setdefault(geneid, [])
		generows.append([ geneid, geneid, exonid, coord, "unknown" ])
	generows = genenametorows.setdefault(genename, [])
	generows.append([ genename, genename, exonid, coord, "unknown" ])

# now that we have the exon ids from the two lists, build a master list based on the ids from IAN_GENES.txt

newlist = csv.writer(open("newexonlist.csv", "w"), delimiter="\t")
seen = Set()
for row in ianreader:
	name, chrom, start, end, strand, geneid = row
	if geneid in seen:
		print "seen gene:", geneid, "before!"
		continue
	if name in seen:
		print "seen gene:", name, "before!"
		continue
	if geneid in genetorows:
		generows = genetorows[geneid]
		seen.add(geneid)
	else:
		if name in genenametorows:
			generows = genenametorows[name]
			seen.add(name)
		else:
			print "could not find", name
			continue
	coord = "chr%s:%s-%s" % (chrom, start, end)
	for row in generows:
		row.insert(2, coord)
		newlist.writerow(row)
