# takes a fasta file containing a genome and a site sequence (like "GTAC")
# makes a bed file representing the location of the sites.

# if you want to just count the sites, use countsite.py

import csv,sys
from Bio import SeqIO


# specific to the ensembl fasta format. They look like this:
# chromosome:NCBI36:10:1:135374737:1
def getChrom(seqid):
    return "chr" + seqid.split(":")[2]

# for the NCBI format. they look like this:
# >gi|65493515|ref|NC_000070.2|NC_000070 Mus musculus chromosome 4
def getChromAlt(seqid):
    fastaId = seqid.split(" ")[-1]
    if fastaId.find("chr") != -1: return fastaId 
    else: return "chr" + fastaId


def dumpSites(outputFilename, genome, site):
    inhandle = SeqIO.parse(open(genome, "rU"), "fasta")
    sitelen = len(site)
    site = site.upper()
    print "Writing bedfile of site locations to file: %s" % (outputFilename)
    writer = csv.writer(open(outputFilename, "w"), delimiter="\t")
    for seqrec in inhandle:
        chromname = getChromAlt(seqrec.description)
        #print "working at chrom", chromname
        seq = seqrec.seq.upper()
        start = 0
        while seq.find(site, start) != -1:
            index = seq.find(site, start)
            writer.writerow([chromname, str(index+1), str(index+sitelen)])
            start = index + 1
    return outputFilename

def getSites(genome, site):
    fastaFile = open(genome, "rU")
    inhandle = SeqIO.parse(fastaFile, "fasta")
    sitelen = len(site)
    site = site.upper()
    sites = []
    siteDict = {}
    positionCounter = 0
    for seqrec in inhandle:
        chromname = getChromAlt(seqrec.description)
        #print "working at chrom", chromnamesiteDict
        seq = seqrec.seq.upper()
        start = 0
        initialPosition = positionCounter
        while seq.find(site, start) != -1:
            index = seq.find(site, start)
            sites.append([chromname, index+1, index+sitelen])
            start = index + 1
            positionCounter += 1
        siteDict[chromname] = [initialPosition, positionCounter]
    fastaFile.flush()
    fastaFile.close()
    return sites, siteDict

def getSitesDict(genome, site):
    fastaFile = open(genome, "rU")
    inhandle = SeqIO.parse(fastaFile, "fasta")
    sitelen = len(site)
    site = site.upper()
    sites = {}
    for seqrec in inhandle:
        chromname = getChromAlt(seqrec.description)
        #print "working at chrom", chromname
        seq = seqrec.seq.upper()
        start = 0
        while seq.find(site, start) != -1:
            index = seq.find(site, start)
            sites[chromname+"_"+str(index+1)] = [chromname, index+1, index+sitelen]
            start = index + 1
    fastaFile.flush()
    fastaFile.close()
    return sites

def obtainSitesBedFile(outputFilename, genome, site):
    return dumpSites(outputFilename, genome, site)

def obtainSitesArray(genome, site):
    return getSites(outputFilename, genome, site)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print len(sys.argv)
        print sys.argv
        print "usage: sitetobed.py <outputfilename> <genome.fa> <siteseq>"
        sys.exit(1)
    outputFilename = sys.argv[1]
    genome = sys.argv[2]
    site = sys.argv[3]
    print outputFilename, genome, site
    dumpSites(outputFilename, genome, site)