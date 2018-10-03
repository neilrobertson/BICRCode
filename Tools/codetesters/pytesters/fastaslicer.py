#!/usr/bin/env python

from Bio import SeqIO, SeqRecord
import sys

def seqSlice(inseqrecord, regstart, globstart, globend):
	# the sequence bit
	baseseq = inseqrecord.seq
	start = globstart - regstart
	end = globend - regstart
	thisslice = baseseq[start:end]
	# the name bit
	newnamelist = inseqrecord.name.split(":")[:3]
	newnamelist.extend([ globstart, globend ])
	newnamelist = [ str(x) for x in newnamelist ]
	newname = ":".join(newnamelist)
	# the new sec record
	seqrec = SeqRecord.SeqRecord(thisslice, newname, "", "")
	return seqrec


fastafile = "tester.fasta"

handle = open(fastafile)
records = [ record for record in SeqIO.parse(handle, "fasta") ]

assert len(records) == 1

record = records[0]

slicestart = 110517481
sliceend = 110517581

seqrecslice = seqSlice(record, slicestart, slicestart, sliceend)


print seqrecslice
