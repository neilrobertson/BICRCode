'''
Created on 28 Feb 2011

@author: mcbryan
'''
import sys
import getopt
from bed.treatment import BedFile
from genemapping import Ensembl
from genemapping.chrmEnds import ChromosomeEnds

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", [
                                                      # command args go here
                                                      "bedfile="
                                                      ])
        
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        print "Usage: main.py  [Space seperated list of gene id sets]"
        sys.exit(2)
    
    
    for o, a in opts:
        if (o=="--bedfile"):
            bedfile = BedFile(a)
    
    intervalTree = bedfile.buildIntervalTree()
    
    chromosomeEnds = ChromosomeEnds("hg18")
    
    genedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")
    
    for gene in genedata:
        # look up data
        if genedata[gene].chr in chromosomeEnds:
            counts = len(intervalTree.getValuesOfOverlappingIntervals(genedata[gene].chr,
                                                                      genedata[gene].start,
                                                                      genedata[gene].end))
            length = genedata[gene].end - genedata[gene].start
            
            countsPerBP = float(counts)/float(length)
            
            # output gene + counts per bp
            print gene,str(countsPerBP)