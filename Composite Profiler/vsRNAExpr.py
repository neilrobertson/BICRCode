'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
from csvfile.indexedcsv import IndexedCSV
from genemapping.Ensembl import EnsemblGenes
import csv
import sequence.genome
from csvfile.genelist import GeneList
import re

from Data import Bed,BedWithControl


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["aSample=","bSample=",
                                                      "aCtrl=","bCtrl=",
                                                      "gene-expression-file=",
                                                      "fccol=",
                                                      "exprcols=",
                                                      "outputfile=",
                                                      "promotorsize=",
                                                      "genelist=",
                                                      "exactmatch",
                                                      "assembly="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    ensemblidcol = "ensemblid"
    fccol = "log2(fold_change)"
    exprcols = ["value_1","value_2"]
    
    extend = 150
    
    exprCSV= None
    outputfile = None
    
    upstreamPromotor = 5000
    downstreamPromotor = 1000
    
    genelists = []
    exactmatch = False
    
    genedata = EnsemblGenes(assembly="hg18")
    
    for opt, value in opts:
        if opt=="--gene-expression-file":
            exprCSV = IndexedCSV(value,keyPos=1)
            print "Expr File:",value
        elif opt=="--fccol":
            fccol = value
        elif opt=="--exprcols":
            exprcols = value.split(",")
        elif opt=="--outputfile":
            outputfile = value
        elif opt=="--promotorsize":
            upstreamPromotor = int(value)
            downstreamPromotor = int(value)
        elif opt=="--aSample":
            aSample = Bed(value,extend=extend)
            print "A Sample:",value
        elif opt=="--bSample":
            bSample = Bed(value,extend=extend)
            print "B Sample:",value
        elif opt=="--genelist":
            genelists.append(GeneList(value))
        elif opt=="--exactmatch":
            exactmatch = True
        elif opt=="--assembly":
            genedata = EnsemblGenes(assembly=value)
            print "Assembly:",value


    # stuff that needs to be done after primary data is loaded
    for o,value in opts:
        if o == "--aCtrl":
            aSample = BedWithControl(value,aSample,extend=extend)
            print "A Ctrl:",value
        if o == "--bCtrl":
            bSample = BedWithControl(value,bSample,extend=extend)
            print "B Ctrl:",value

    def getRegion(a,chrm,start,stop):
        l = a.valuesBehaviour(chrm,a.getValues(chrm, start, stop),start,stop)
        assert len(l) == 1
        return l[0]

    def regionDifference(a,b,chrm,start,stop):
        aValue = getRegion(a,chrm,start,stop)
        bValue = getRegion(b,chrm,start,stop)
        return bValue - aValue

    assert aSample != None
    assert bSample != None

    assert exprCSV != None
    assert outputfile != None       
  
    outputcsv = csv.writer(open(outputfile,"w"),delimiter="\t")
    
    headerrow = ["ensemblid","genename","logfc"]
    headerrow.extend(exprcols)
    headerrow.extend(["promotorDiff","regionA","regionB"])
    headerrow.extend(["genebodyDiff","genebodyA","genebodyB"])
    headerrow.extend(["promoterGenebodyDiff","promoterGenebodyA","promoterGenebodyB"])
    
    for genelist in genelists:
        headerrow.append(genelist.getFriendlyName())
    outputcsv.writerow(headerrow)
    
    # make a list of genenames in each genelist and store in the genelist
    for genelist in genelists:
        genelist.genenames = set()
        for e in genelist:
            if e in genedata:
                genelist.genenames.add(genedata[e].name)
    
    for testid in exprCSV:
        try:
            ensembl = exprCSV[testid][ensemblidcol]
            
            if ensembl not in genedata:
                # RNA-Seq is not a single unique gene, skip
                continue
        
            logFC = exprCSV[testid][fccol]
                    
            promotorstart, promotorstop = genedata[ensembl].getPromotorRegion(upstreamPadding = upstreamPromotor, downstreamPadding = downstreamPromotor)          
    
            row = [ensembl,genedata[ensembl].name, logFC]
            
            exprs = []
            for exprcol in exprcols:
                exprs.append(exprCSV[testid][exprcol]) 
            row.extend(exprs)
            
            promotorDiff = regionDifference(aSample,bSample,genedata[ensembl].chr,promotorstart,promotorstop)
            row.append(promotorDiff)
            
            promotorA = getRegion(aSample,genedata[ensembl].chr,promotorstart,promotorstop)
            promotorB = getRegion(bSample,genedata[ensembl].chr,promotorstart,promotorstop)
            
            row.extend([promotorA,promotorB])
            
            
            genebodyDiff = regionDifference(aSample,bSample,genedata[ensembl].chr,genedata[ensembl].start,genedata[ensembl].end)
            row.append(genebodyDiff)
            
            genebodyA = getRegion(aSample,genedata[ensembl].chr,genedata[ensembl].start,genedata[ensembl].end)
            genebodyB = getRegion(bSample,genedata[ensembl].chr,genedata[ensembl].start,genedata[ensembl].end)

            row.extend([genebodyA,genebodyB])
            
            
            promoterGenebodyDiff = regionDifference(aSample,bSample,genedata[ensembl].chr,genedata[ensembl].start,genedata[ensembl].end)
            row.append(genebodyDiff)
            
            if promotorstart < genedata[ensembl].start:
                promoterGenebodyDiff = regionDifference(aSample,bSample,genedata[ensembl].chr,promotorstart,genedata[ensembl].end)
                row.append(genebodyDiff)
                promoterGenebodyA = getRegion(aSample,genedata[ensembl].chr,promotorstart,genedata[ensembl].end)
                promoterGenebodyB = getRegion(bSample,genedata[ensembl].chr,promotorstart,genedata[ensembl].end)
                row.extend([promoterGenebodyA,promoterGenebodyB])

            elif promotorstop > genedata[ensembl].end:
                promoterGenebodyDiff = regionDifference(aSample,bSample,genedata[ensembl].chr,genedata[ensembl].start,promotorstop)
                row.append(genebodyDiff)
                promoterGenebodyA = getRegion(aSample,genedata[ensembl].chr,genedata[ensembl].start,promotorstop)
                promoterGenebodyB = getRegion(bSample,genedata[ensembl].chr,genedata[ensembl].start,promotorstop)
                row.extend([promoterGenebodyA,promoterGenebodyB])

            
            for genelist in genelists:
                                
                matchespattern = (ensembl in genelist.seengenes) or (genedata[ensembl].name in genelist.genenames)
                
                if not exactmatch:
                    for pattern in genelist:
                        if matchespattern:
                            break
                        if re.match(pattern,genedata[ensembl].name):
                            matchespattern = True
                        
                row.append("1" if matchespattern else "0")
            
            outputcsv.writerow(row)
        except sequence.genome.UnknownChromosomeException:
            continue