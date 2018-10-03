'''
Created on 18 Mar 2011

@author: mcbryan
'''

import sys
import getopt
import csv
import collections

from genemapping import Ensembl

if __name__ == '__main__':
    
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "G:o:", [
                                                      # command args go here
                                                      "gtf=",
                                                      "output="
                                                      ])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        print "Usage: main.py  [Space seperated list of gene id sets]"
        sys.exit(2)
        
    for o, a in opts:
        if (o=="-G") or (o=="--gtf"):
            gtfFile = azz
        elif (o=="-o") or (o=="--output"):
            outputFile = a
        else:
            print "Unknown parameter: "+o+" "+a
            sys.exit(2)
    
    genedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")
    
    gtfReader = csv.reader(open(gtfFile,"r"),delimiter="\t")

    cuffGeneIdsToEnsemblGeneIds = collections.defaultdict(set)
    
    missing= set()
    for row in gtfReader:
        detailsColumn = row[8]
        details = dict(item.replace("\"","").split(" ") for item in detailsColumn.split("; "))
                
        if "nearest_ref" not in details:
            continue
        
        if details["nearest_ref"] not in genedata.transcripts:
            print "Missing ref:"+details["nearest_ref"]
            missing.add(details["nearest_ref"])
            continue
        
        geneid = details["gene_id"]
        nearest_ref = details["nearest_ref"]
        
        cuffGeneIdsToEnsemblGeneIds[geneid].add(genedata.transcripts[nearest_ref].geneid)
    
    multipletotal = 0
    total=0
    for cuffGeneId in cuffGeneIdsToEnsemblGeneIds:
        total+=1
        if len(cuffGeneIdsToEnsemblGeneIds[cuffGeneId]) > 1:
            multipletotal+=1
    
    print "Multiple:"+str(multipletotal)+"/"+str(total)
    
    output = open(outputFile,"w")
    
    gtfReader = csv.reader(open(gtfFile,"r"),delimiter="\t",)
    
    for row in gtfReader:
        detailsColumn = row[8]
        details = dict(item.replace("\"","").split(" ") for item in detailsColumn.split("; "))
        
        geneid = details["gene_id"]
        
        if len(cuffGeneIdsToEnsemblGeneIds[geneid])==1:
            row[8] = str(row[8].replace(geneid,list(cuffGeneIdsToEnsemblGeneIds[geneid])[0]))
            #print row[8]
        print >> output, "\t".join(row)
        