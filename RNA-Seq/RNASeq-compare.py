'''
Created on 5 Aug 2010

@author: mcbryan
'''



import getopt
import sys
import csv

from csvfile.indexedcsv import IndexedCSV

from genemapping import UCSC
from genemapping import Ensembl
from affy.NetAffxAnnotation import NetAffxAnnotation

from matplotlib.pyplot import *


if __name__ == '__main__':
    
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["microarray=","rnaseq=", "output="])
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
    finalentries = 0
    
    with open(outputFileLoc,"w") as outputfile:
        outcsv = csv.writer(outputfile,delimiter="\t")
    
        for affy in affyfile:
            
            allaffys+=1
            
            affypval = affyfile[affy]["tt"]
        
            if float(affypval)>0.05:
                continue
            
            ensemblGenes = affyannotation.getValues(affy,"Ensembl")
            
            if len(ensemblGenes) != 1:
                continue # we only want unique probes here
            else:
                if ensemblGenes[0] == "---":
                    continue
                ensemblGene = ensemblgenedata[ensemblGenes[0]]
            
            affy_fc = affyfile[affy]["log2(fc)"]   
            
            ucscTranscripts = set()
            for ensemblTranscript in ensemblGene:
                ucscTranscripts.update(ucscgenedata.getTranscriptsForEnsembl(ensemblTranscript))
            
            clusters = []
            for ucscTranscript in ucscTranscripts:
                clusters.append(ucscgenedata[ucscTranscript].clusterid)
            
            #print ensemblGene,affy_fc, ucscTranscripts, clusters
            
            if len(clusters) != 1:
                continue # we only want unique mappings here
            else:
                clusterid = clusters[0]
            
            if clusterid not in rnaseqfile:
                continue
            
            baseMean = float(rnaseqfile[clusterid]["baseMean"])
            
            rna_fc = rnaseqfile[clusterid]["log2FoldChange"]
            
            if affy_fc == "NA" or rna_fc in ["NA","Inf","-Inf"] or baseMean < 20.0:
                continue
        
            rnapval = rnaseqfile[clusterid]["pval"]
        
            if float(rnapval)>0.05:
                continue
            
            print [affy_fc,rna_fc]
            
            outcsv.writerow([affy_fc,rna_fc])
            
            plotx.append(float(affy_fc))
            ploty.append(float(rna_fc))
            
            finalentries += 1
    
    print allaffys
    print finalentries
    
    plot(plotx, ploty, 'o')
    show()