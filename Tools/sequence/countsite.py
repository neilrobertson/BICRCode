from Bio import SeqIO

# takes a fasta file containing a whole genome and a sequence (site)
# counts the number of sites

def chromkey(x):
	try:
		return int(x[3:])
	except ValueError:
		return x[3:]

# specific to the ensembl fasta format. They look like this:
# chromosome:NCBI36:10:1:135374737:1
def getChrom(seqid):
	return "chr" + seqid.split(":")[2]

def countSites(genome, site):
	sitecount = {}
	inhandle = SeqIO.parse(open(genome, "rU"), "fasta")
	for seqrec in inhandle:
		chromname = getChrom(seqrec.id)
		print "working at chrom", chromname
		count = seqrec.seq.count(site)
		sitecount[chromname] = count
	return sitecount

if __name__ == "__main__":
	import sys
	if len(sys.argv) != 3:
		print "usage: countsite.py <genome.fasta> <sitesequence>"
		sys.exit(1)
	genome = sys.argv[1]
	site = sys.argv[2]
	sitecount = countSites(genome, site)
	chroms = sitecount.keys()
	chroms.sort(key=chromkey)
	for chrom in chroms:
		print "%s: %d" % (chrom, sitecount[chrom])
	print "=="
	print "total: %d" % (sum(sitecount.values()))
