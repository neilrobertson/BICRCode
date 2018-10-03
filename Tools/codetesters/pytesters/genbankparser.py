#!/usr/bin/env python

from Bio import SeqIO
import csv

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

gbfile = "/home/pzs/primerdesign/promoters/skurukutipromo.gb"
#gbfile = "tester.gb"
outfile = "promoterseqs.csv"

records = SeqIO.parse(open(gbfile, "rb"), "genbank")

writer = csv.writer(open(outfile, "wb"))

for record in records:
	descmap = {}
	descs = record.description.split("|")
	for desc in descs:
		try:
			key, val = desc.split("=")
			descmap[key.strip()] = val.strip()
		except ValueError:
			print "problem with description item", desc
			print "in record", record
	name = descmap["acc"]
	coord = "chr%s:%s-%s" % (descmap["chr"], descmap["start"], descmap["end"])
	strand = descmap["str"]
	if strand == "(-)":
		row = [ name, coord, record.seq.reverse_complement() ]
	elif strand == "(+)":
		row = [ name, coord, record.seq ]
	else:
		print "unknown strand type:", strand
		print "in record", record
		continue
	writer.writerow(row)

	
