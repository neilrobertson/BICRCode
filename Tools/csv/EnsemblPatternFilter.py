'''
Created on 11 Feb 2011

@author: mcbryan
'''
from genemapping.chrmEnds import ChromosomeEnds
from genemapping import Ensembl
import sys
import getopt
import re
from csvfile.genelist import GeneList

unique = True


if __name__ == '__main__':
    
    patterns = None
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "m:", [])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
    
    for opt, value in opts:
        if opt=="-m":
            patterns = GeneList(value)
            
    assert patterns != None
    
    chromosomeEnds = ChromosomeEnds("hg18")
    
    genedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")
    
    for file in args:
        currentfile = GeneList(file)
        
        for gene in currentfile:
            name = genedata[gene].name
            match = False
            for pattern in patterns:
                if re.match(pattern,genedata[gene].name):
                    match = True
                    break
            if match:
                print "Y",gene,name
            else:
                print "N",gene,name