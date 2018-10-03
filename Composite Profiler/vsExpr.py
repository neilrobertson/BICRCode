'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
from affy.NetAffxAnnotation import NetAffxAnnotation
from csvfile.indexedcsv import IndexedCSV
import collections
import math
from genemapping.Ensembl import EnsemblGenes
import csv
import sequence.genome
from sequence.genome import Genome
from csvfile.genelist import GeneList
import re
from datastructures.memoized import memoized

from Data import Bed,BedWithControl


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["aSample=","bSample=",
                                                      "aCtrl=","bCtrl=",
                                                      "affy-gene-expression=",
                                                      "affyfccol=",
                                                      "affyexprcol=",
                                                      "outputfile=",
                                                      "promotorsize=",
                                                      "genelist=",
                                                      "build=",
                                                      "exactmatch"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    affyfile= None
    outputfile = None
    affyfccol = None
    affyexprcol = None
    upstreamPromotor = 2000
    downstreamPromotor = 2000
    exactmatch = False
    
    genelists = []
    
    for opt, value in opts:
        if opt=="--affy-gene-expression":
            affyfile = value
        elif opt=="--affyfccol":
            affyfccol = value
        elif opt=="--affyexprcol":
            if affyexprcol == None:
                affyexprcol = [value]
            else:
                affyexprcol.append(value)
        elif opt=="--outputfile":
            outputfile = value
        elif opt=="--promotorsize":
            upstreamPromotor = int(value)
            downstreamPromotor = int(value)
        elif opt=="--aSample":
            aSample = Bed(value)            
        elif opt=="--bSample":
            bSample = Bed(value)
        elif opt=="--genelist":
            genelists.append(GeneList(value))
        elif opt=="--build":
            build = value
        elif opt=="--exactmatch":
            exactmatch = True

    # stuff that needs to be done after primary data is loaded
    for o,value in opts:
        if o == "--aCtrl":
            aSample = BedWithControl(value,aSample)
        if o == "--bCtrl":
            bSample = BedWithControl(value,bSample)

    def getAvg(values):
        if len(values)==0:
            return None
        return math.fsum(values) / float(len(values))

    @memoized
    def getRegion(a,chrm,start,stop):
        return getAvg(a.valuesBehaviour(chrm,a.getValues(chrm, start, stop),start,stop))

    def regionDifference(a,b,chrm,start,stop):
        aValue = getRegion(a,chrm,start,stop)
        bValue = getRegion(b,chrm,start,stop)
        return bValue - aValue

    assert aSample != None
    assert bSample != None

    assert affyfile != None
    assert affyfccol != None
    assert affyexprcol != None
    assert outputfile != None
    
    print affyexprcol
    
    genedata = EnsemblGenes(assembly=build)
    
    genome = Genome(genomeBuild = build)
    
    affyannotation = NetAffxAnnotation(genome = build, cdfname="HG-U133_Plus_2")
    
    affyCSV = IndexedCSV(affyfile)
    affyEnsemblLogFCs = collections.defaultdict(list)
    affyEnsemblExprs = collections.defaultdict(list)
    
    # make a list of genenames in each genelist and store in the genelist
    for genelist in genelists:
        genelist.genenames = set()
        for e in genelist:
            if e in genedata:
                genelist.genenames.add(genedata[e].name)
    
    for affy in affyCSV:
        ensembls = affyannotation.getValues(affy, "Ensembl")
        if len(ensembls)==1:
            affyFC = float(affyCSV[affy][affyfccol])
            affylogFC = math.log(affyFC) if affyFC > 0.0 else math.log(abs(affyFC))*-1.0
            affyEnsemblLogFCs[ensembls[0]].append(affylogFC)
            
            affyexpr = {}
            for col in affyexprcol:
                affyexpr[col]=float(affyCSV[affy][col])
            
            affyEnsemblExprs[ensembls[0]].append(affyexpr)
        
  
    outputcsv = csv.writer(open(outputfile,"w"),delimiter="\t")
    
    headerrow = ["ensemblid","genename","affylogfc"]
    headerrow.extend(affyexprcol)
    headerrow.extend(["promotorDiff","regionA","regionB"])
    for genelist in genelists:
        headerrow.append(genelist.getFriendlyName())    
    outputcsv.writerow(headerrow)
    
    for ensembl in affyEnsemblLogFCs:
        try:
            affyavglogFC = math.fsum(affyEnsemblLogFCs[ensembl])/float(len(affyEnsemblLogFCs[ensembl]))
            
            affyavgexpr = collections.defaultdict(list)
            for affyexpr in affyEnsemblExprs[ensembl]:
                for col in affyexprcol:
                    affyavgexpr[col].append(affyexpr[col])
            
            promotorstart, promotorstop = genedata[ensembl].getPromotorRegion(upstreamPadding = upstreamPromotor, downstreamPadding = downstreamPromotor)          
            
            
            row = [ensembl,genedata[ensembl].name, affyavglogFC]
            
            for col in affyexprcol:
                row.append(math.fsum(affyavgexpr[col])/float(len(affyavgexpr[col])))
                
                
            promotorDiff = regionDifference(aSample,bSample,genedata[ensembl].chr,promotorstart,promotorstop)                
            row.append(promotorDiff)
            
            regionA = getRegion(aSample,genedata[ensembl].chr,promotorstart,promotorstop)
            regionB = getRegion(bSample,genedata[ensembl].chr,promotorstart,promotorstop)
            
            row.extend([regionA,regionB])
            
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
        except KeyError:
            continue