'''
Created on 11 Feb 2011

@author: mcbryan
'''
from genemapping import Ensembl
import sys
import getopt
from csvfile.genelist import GeneList

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", [])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
    genedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")
    
    for file in args:
        currentfile = GeneList(file)
        
        for gene in currentfile:
            print str(genedata[gene].chr) + "\t" + str(genedata[gene].start) + "\t" + str(genedata[gene].end)