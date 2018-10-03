#!/usr/bin/env python

import sys

# lookup sequence by coordinates.
# this first implementation uses fasta files

class GeneSeqLookup:
	def __init__(self):
		self._chromseqs = {}
		
	def stripPath(self, pathname):
		from os.path import basename
		return basename(pathname).split(".")[0]

	def initialiseFromFastaDirectory(self, dirname):
		from glob import glob
		from Bio import SeqIO
		print >> sys.stderr, "initialising gene sequence lookup from directory", dirname
		fastafiles = glob(dirname + "/*.fasta")
		if len(fastafiles) == 0:
			fastafiles = glob(dirname + "/*.fa")
		for ffile in fastafiles:
			chrname = self.stripPath(ffile)
			chrname = chrname.replace("chr", "")
			# skip the whole-genome fasta file
			if chrname.startswith("all"):
				continue
			handle = open(ffile)
			seq_records = [ sr for sr in SeqIO.parse(handle, "fasta") ]
			assert len(seq_records) == 1
			self._chromseqs[chrname] = seq_records[0].seq
			
	def getSeq(self, chrom, start, end):
		chrom = chrom.replace("chr", "")
		chromseq = self._chromseqs[chrom]
		return chromseq[start:end]

	def getGCProportion(self, chrom, start, end):
		seq = self.getSeq(chrom, start, end)
		gccontent = seq.count("G") + seq.count("C")
		return float(gccontent) / float(len(seq))


if __name__ == "__main__":
	gsl = GeneSeqLookup()
	gendir = "/home/pzs/genebuilds/mouse/"
	gsl.initialiseFromFastaDirectory(gendir)
	print gsl.getSeq("1", 3, 6)
	print gsl.getGCProportion("1", 3, 6)
	print gsl.getSeq("2", 3, 10)
	print gsl.getGCProportion("2", 3, 10)
	print gsl.getSeq("5", 54434082, 54434106)
	print gsl.getSeq("12", 17411777, 17411856)
