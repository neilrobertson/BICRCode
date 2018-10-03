# takes a fasta file containing a genome and a site sequence (like "GTAC")
# makes a bed file representing the location of the sites.

# if you want to just count the sites, use countsite.py

import sys
from Bio import SeqIO

# specific to the ensembl fasta format. They look like this:
# chromosome:NCBI36:10:1:135374737:1
def getChrom(seqid):
	return "chr" + seqid.split(":")[2]

def dumpSites(genome, site):
	inhandle = SeqIO.parse(open(genome, "rU"), "fasta")
	sitelen = len(site)
	site = site.upper()
	for seqrec in inhandle:
		chromname = getChrom(seqrec.id)
		#print "working at chrom", chromname
		seq = seqrec.seq
		start = 0
		while seq.find(site, start) != -1:
			index = seq.find(site, start)
			outlist = [ chromname, str(index+1), str(index+sitelen) ]
			print "\t".join(outlist)
			start = index + 1

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: sitetobed.py <genome.fa> <siteseq>"
		sys.exit(1)
	genome = sys.argv[1]
	site = sys.argv[2]
	dumpSites(genome, site)
