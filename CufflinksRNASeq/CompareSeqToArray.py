'''
Created on 21 Mar 2011

@author: mcbryan
'''

import sys
import getopt
import math
import matplotlib
import matplotlib.pyplot as plt
from csvfile.indexedcsv import IndexedCSV
from affy.NetAffxAnnotation import NetAffxAnnotation
import collections

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "r:a:", ["affypvalue=","affyfc=","rnalogfc=","rnasignificant=","rnakey="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    affypvalue = "BH-fdr"
    affyfc = "fc"
    
    rnalogfc = "ln(fold_change)"
    rnasig = "significant"
    
    rnakeypos = 0
    
    for o, a in opts:
        if o=="-r":
            rnaseqfile = a
        elif o=="-a":
            affyfile = a
        elif o=="--affypvalue":
            affypvalue = a
        elif o=="--affyfc":
            affyfc = a
        elif o=="--rnalogfc":
            rnalogfc = a
        elif o=="--rnasignificant":
            rnasig = a
        elif o=="--rnakey":
            rnakeypos = int(a)
    
    # read in gene_exp.diff as indexedcsv
    
    affyannotation = NetAffxAnnotation(genome = "hg18", array="HG-U133_Plus_2", version="29")
    
    affyCSV = IndexedCSV(affyfile)
    rnaseqCSV = IndexedCSV(rnaseqfile,keyPos=rnakeypos)
    
    fig = plt.figure(figsize=(12, 12), dpi=100)
    ax = fig.add_subplot(111)
    
    affyEnsemblLogFCs = collections.defaultdict(list)
    affyEnsemblSig = collections.defaultdict(bool)
    
    for affy in affyCSV:
        ensembls = affyannotation.getValues(affy, "Ensembl")
        if len(ensembls)==1:
            affyFC = float(affyCSV[affy][affyfc])
            affylogFC = math.log(affyFC) if affyFC > 0.0 else math.log(abs(affyFC))*-1.0
            affyEnsemblLogFCs[ensembls[0]].append(affylogFC)
            
            if float(affyCSV[affy][affypvalue])<0.05:
                affyEnsemblSig[ensembls[0]] = True
            else:
                affyEnsemblSig[ensembls[0]] = False
    
    x = []
    y = []
    cols = []
    
    both = 0
    affyonly = 0
    rnaonly = 0
    neither = 0
    
    for ensembl in affyEnsemblLogFCs:

        if ensembl in rnaseqCSV:
            # store it
            
            
            rnaseqlogFC = float(rnaseqCSV[ensembl][rnalogfc])
            affylogFC = math.fsum(affyEnsemblLogFCs[ensembl])/float(len(affyEnsemblLogFCs[ensembl]))

            if rnaseqlogFC > 5.0:
                rnaseqlogFC = 5.0                
            elif rnaseqlogFC < -5.0:
                rnaseqlogFC = -5.0
            
#            x.append(rnaseqlogFC)
#            y.append(affylogFC)
            
            if affyEnsemblSig[ensembl]==True and rnaseqCSV[ensembl][rnasig]=="yes":
                both +=1
                x.append(rnaseqlogFC)
                y.append(affylogFC)
                cols.append("b")
            elif affyEnsemblSig[ensembl]==True:
                affyonly +=1
                x.append(rnaseqlogFC)
                y.append(affylogFC)
                cols.append("r")
            elif rnaseqCSV[ensembl][rnasig]=="yes":
                rnaonly+=1
                x.append(rnaseqlogFC)
                y.append(affylogFC)
                cols.append("g")
            else:
                neither+=1
                x.append(rnaseqlogFC)
                y.append(affylogFC)
                cols.append("w")
    
    print "Both:"+str(both)
    print "Affy only:"+str(affyonly)
    print "RNA only:"+str(rnaonly)
    print "Neither:"+str(neither)
    
    ax.scatter(x, y, c=cols, marker='o')
    plt.xlabel("RNA-seq FC (log)")
    plt.ylabel("Affy FC (log)")
    plt.xlim([-6,6])
    plt.ylim([-2,2])
    plt.show()