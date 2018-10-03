#!/usr/bin/env python

import sys
import getopt
from genemapping import Ensembl
from csvfile.genelist import GeneList
import re

try:
    opts, args = getopt.getopt(sys.argv[1:], "g:", [])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    sys.exit(2)

UPSTREAM_PROMOTOR_DIST = 5000

genelists = []

for o, a in opts:
    if o=="-g":
        genelists.append(GeneList(a))

assert len(genelists) > 0

genedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")

for genelist in genelists:
    print genelist.getFullName()
    for gene in genedata:
        # does the gene match the pattern
        for pattern in genelist:
            if re.match(pattern,genedata[gene].name):
                start,stop = genedata[gene].getGeneWithPromotor(upstreamPadding = UPSTREAM_PROMOTOR_DIST)
                
                print genedata[gene].id, genedata[gene].name, genedata[gene].chr, start, stop
                break