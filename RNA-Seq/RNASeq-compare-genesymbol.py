'''
Created on 5 Aug 2010

@author: mcbryan
'''



import getopt
import sys
import math
import csv

from csvfile.indexedcsv import IndexedCSV

from genemapping import UCSC
from genemapping import Ensembl
from affy.NetAffxAnnotation import NetAffxAnnotation

from matplotlib.pyplot import *

import collections

if __name__ == '__main__':
    
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["microarray=","rnaseq=","output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    affyfilelocation = None
    rnaseqfilelocation = None
    outputFileLoc = None
    
    for o, a in opts:
        if (o=="--microarray"):
            affyfilelocation = a
        elif (o=="--rnaseq"):
            rnaseqfilelocation = a
        elif (o=="--output"):
            outputFileLoc = a
    
    assert affyfilelocation != None and rnaseqfilelocation != None
    assert outputFileLoc != None
    
    ucscgenedata = UCSC.UCSCTranscripts(assembly="hg18")
    ensemblgenedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")
    
    affyannotation = NetAffxAnnotation()
    
    affyfile = IndexedCSV(affyfilelocation)
    rnaseqfile = IndexedCSV(rnaseqfilelocation,keyPos=1)
    
    print "Read Files"
    
    plotx = []
    ploty = []
    
    allaffys = 0
    allrnaclusters = 0
    finalentries = 0
    
    
    
    affySymbols = collections.defaultdict(list)
    rnaSymbols = collections.defaultdict(list)
    
    for affy in affyfile:
        
        allaffys+=1
        
        symbol = affyfile[affy]["Gene symbol"]
        
        affy_fc = affyfile[affy]["log2(fc)"]
        #affy_fc = affyfile[affy]["fc"]
        
        pval = affyfile[affy]["tt"]
        
        if float(pval)>0.05:
            continue
        
        affySymbols[symbol].append(float(affy_fc))
        

        
    for clusterid in rnaseqfile:
        
        allrnaclusters += 1
        
        symbol = rnaseqfile[clusterid]["Symbol"]
        
        rna_fc = rnaseqfile[clusterid]["log2FoldChange"]
        #rna_fc = rnaseqfile[clusterid]["foldChange"]
        
        baseMean = float(rnaseqfile[clusterid]["baseMean"])
        
        if rna_fc in ["NA","Inf","-Inf"] or baseMean < 20.0:
            continue
        
        pval = rnaseqfile[clusterid]["pval"]
        
        if float(pval)>0.05:
            continue
        
        rnaSymbols[symbol].append(float(rna_fc))
    
    with open(outputFileLoc,"w") as outputfile:
        outcsv = csv.writer(outputfile,delimiter="\t")
        for symbol in affySymbols:
            if symbol in rnaSymbols:
                
                if len(affySymbols[symbol]) > 1 or len(rnaSymbols[symbol]) > 1:
                    continue
                
                affy_fc = affySymbols[symbol][0]
                rna_fc = rnaSymbols[symbol][0]
                
                #affy_fc = math.fsum(affySymbols[symbol])/float(len(affySymbols[symbol]))
                #rna_fc = math.fsum(rnaSymbols[symbol])/float(len(rnaSymbols[symbol]))
            
                plotx.append(float(affy_fc))
                ploty.append(float(rna_fc))
                
                outcsv.writerow([affy_fc,rna_fc])
            
                finalentries += 1
    
    print len(affySymbols)
    print len(rnaSymbols)
    print finalentries
    
    plot(plotx, ploty, 'o')
    show()