'''
Created on 11 Feb 2011

@author: mcbryan
'''
from genemapping.chrmEnds import ChromosomeEnds
from genemapping import Ensembl
import fileinput

unique = True

if __name__ == '__main__':
    chromosomeEnds = ChromosomeEnds("hg18")
    
    genedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")

    def valid(gene):
        if getChr(gene) not in chromosomeEnds:
            return False
        if gene not in genedata:
            print gene + " not in genedata"
            return False
        return True

    def getChr(gene):
        return genedata[gene].chr
            
    seengenes = set()

    for line in fileinput.input():
        geneid = line.strip()
        if valid(geneid):
            if not unique or geneid not in seengenes:
                print geneid
            seengenes.add(geneid)