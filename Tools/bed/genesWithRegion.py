#!/usr/bin/env python

# prints out only genes (with promotor) that are overlapped by a region


import os
import sys

# add executing directory as path
sys.path.append(sys.path[0])
sys.path.append(os.path.expanduser("~/mount/repository/shared/baselib"))

import csv
import getopt
import numpy
import math
from genemapping import Ensembl
from bed.treatment import Bed as BedIntervalTree
from csvfile.genelist import GeneList

# treatment behaviours

if __name__ == "__main__":
    
    genedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:f:", [])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    # add executing directory as part of path
    sys.path.append(sys.path[0])

    tsspadding = 2000
    filterGeneList = None

    for o, a in opts:
        if o == "-a":
            bTreatment= BedIntervalTree(a)
        if o == "-f":
            filterGeneList = GeneList(a)
    
    countValid = 0
    
    genes = set()
    
    for gene in genedata:
        chr = genedata[gene].chr
        
        if chr not in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                       "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", 
                       "chr21", "chr22", "chrX", "chrY","chrM"]:
            continue
        
        countValid += 1
        
        (start, stop) = genedata[gene].getGeneWithPromotor(upstreamPadding=tsspadding)

        intervals = bTreatment.getOverlappingIntervals(chr, start, stop)
        
        if len(intervals)>0:
            genes.add(gene)
        
    #print countValid
        
    for gene in genes:
        if filterGeneList == None or gene in filterGeneList:
            print gene
