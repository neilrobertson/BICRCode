'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
from datastructures.AggregateTree import AggregateTree
from affy.NetAffxAnnotation import NetAffxAnnotation
from csvfile.indexedcsv import IndexedCSV
import collections
import math
from genemapping.Ensembl import EnsemblGenes
import csv
import os
from bed.treatment import ExtendedBed
import sequence.genome
from sequence.genome import Genome

# https://github.com/brentp/fishers_exact_test 
from fisher import cfisher as fisher_exact

def fisherExact(lmeth,lunmeth,rmeth,runmeth):
    
    p = fisher_exact.pvalue(lmeth,lunmeth,rmeth,runmeth).two_tail
    #return -10.0*math.log10(p)
    return p

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["cpg-meth=","affy-gene-expression=","affyfccol=","affyexprcol=","outputfile=","promotorsize=","affypcol="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    methdata = None
    affyfile= None
    outputfile = None
    
    affyfccol = None
    affyexprcol = None
    affypcol = None
    
    upstreamPromotor = 5000
    downstreamPromotor = 1000
    
    for o, a in opts:
        if o=="--cpg-meth":
            methdatafile = a
        elif o=="--affy-gene-expression":
            affyfile = a
        elif o=="--affyfccol":
            affyfccol = a
        elif o=="--affyexprcol":
            affyexprcol = a
        elif o=="--affypcol":
            affypcol = a
        elif o=="--outputfile":
            outputfile = a
        elif o=="--promotorsize":
            upstreamPromotor = int(a)
            downstreamPromotor = int(a)
    
    assert methdatafile != None
    assert affyfile != None
    assert affyfccol != None
    assert affyexprcol != None
    assert outputfile != None
    
    genedata = EnsemblGenes(assembly="hg18")
    
    genome = Genome(genomeBuild = "hg18")
    
    affyannotation = NetAffxAnnotation(genome = "hg18", cdfname="HG-U133_Plus_2")
    
    cpgIslands = ExtendedBed(os.path.expanduser("~/mount/publicdata/hg18/cpgislands/cpgislands-0-index.bed"))
    
    affyCSV = IndexedCSV(affyfile)
    affyEnsemblLogFCs = collections.defaultdict(list)
    affyEnsemblExprs = collections.defaultdict(list)
    affyEnsemblPvalues = collections.defaultdict(list)
    
    for affy in affyCSV:
        ensembls = affyannotation.getValues(affy, "Ensembl")
        if len(ensembls)==1:
            affyFC = float(affyCSV[affy][affyfccol])
            affylogFC = math.log(affyFC) if affyFC > 0.0 else math.log(abs(affyFC))*-1.0
            affyEnsemblLogFCs[ensembls[0]].append(affylogFC)
            
            affyexpr = float(affyCSV[affy][affyexprcol])
            affyEnsemblExprs[ensembls[0]].append(affyexpr)
            
            affyp = float(affyCSV[affy][affypcol])
            affyEnsemblPvalues[ensembls[0]].append(affyp)

    methdata = AggregateTree(methdatafile)
    
    def methPercentageDiff(values):
        lmethTotal = 0
        lunmethTotal = 0
        rmethTotal = 0
        runmethTotal = 0
        
        for (lmeth,lunmeth,rmeth,runmeth) in values:
            lmethTotal += lmeth
            lunmethTotal += lunmeth
            rmethTotal += rmeth
            runmethTotal += runmeth
            
        lmeth = float(lmethTotal)
        lunmeth = float(lunmethTotal)
        
        # arbitrary cutoff for regions with very few reads
        if lmeth + lunmeth <= 10.0:
            return "", ""
        
        rmeth = float(rmethTotal)
        runmeth = float(runmethTotal)
        
        # arbitrary cutoff for regions with very few reads
        if rmeth + runmeth <= 10.0:
            return "", ""
        
        lpercentage = lmeth / (lmeth+lunmeth)
        rpercentage = rmeth / (rmeth+runmeth)
        
        pvalue = fisherExact(lmethTotal,lunmethTotal,rmethTotal,runmethTotal)
        
        return rpercentage - lpercentage, pvalue
  
  
    def sequenceGCandCpGcontent(sequence):
      
        assert len(sequence) > 0
        
        Ns = sequence.count('N')
        Gs = sequence.count('G')
        Cs = sequence.count('C')
        
        cpgs = sequence.count('CG')
        
        effectiveSeqLen = len(sequence) - Ns
        
        if effectiveSeqLen == 0: # all N's
            gcContent = 0.0
            cpgContent = 0.0
            cpgRatio = 0.0
        else:
            gcContent = float(Gs+Cs) / float(effectiveSeqLen)
            
            cpgContent = float(cpgs) / float(effectiveSeqLen)
            
            if cpgs>0:
                cpgRatio = (float(cpgs) * effectiveSeqLen) / float(Gs*Cs) # from Schubeler
            else:
                cpgRatio = 0.0 # avoid div by 0 error if Gs+Cs == 0, shortcut cpgs == 0 then cpgRatio = 0
    
        return Gs, Cs, cpgs, effectiveSeqLen, gcContent, cpgContent, cpgRatio
  
    
    def classifyCpGContent(sequence):
        Gs, Cs, cpgs, effectiveSeqLen, gcContent, cpgContent, cpgRatio = sequenceGCandCpGcontent(sequence)
        
        hcp = False
        icp = False
        if gcContent >= 0.55 and cpgRatio >= 0.75:
            hcp = True
        elif cpgRatio >= 0.48:
            icp = True
            
        return "high" if hcp else ("intermediate" if icp else "low")
  
  
    outputcsv = csv.writer(open(outputfile,"w"),delimiter="\t")
    outputcsv.writerow(["ensemblid","genename","cpgisland","affylogfc","affyfc","affyavgexpr","affyminp","promotorMethDiff","promotorP","genebodyMethDiff","genebodyP","promotorCpGClassification","genebodyCpGClassification"])
    
    for ensembl in affyEnsemblLogFCs:
        # look up meth value
        
        try:
            affyavglogFC = math.fsum(affyEnsemblLogFCs[ensembl])/float(len(affyEnsemblLogFCs[ensembl]))
            
            affyavgfc = math.exp(affyavglogFC) if affyavglogFC >= 0 else -1*math.exp(abs(affyavglogFC))
            
            affyavgexpr = math.fsum(affyEnsemblExprs[ensembl])/float(len(affyEnsemblExprs[ensembl]))
            affyminp = min(affyEnsemblPvalues[ensembl])
            
            promotorstart, promotorstop = genedata[ensembl].getPromotorRegion(upstreamPadding = upstreamPromotor, downstreamPadding = downstreamPromotor)
            
            incpg = len(cpgIslands.getValuesInRange(genedata[ensembl].chr, promotorstart, promotorstop))!=0            
    
            promotorMethDiff,promotorPvalue = methPercentageDiff(methdata.getValuesInRange(genedata[ensembl].chr,promotorstart,promotorstop))
            geneBodyMethDiff,genebodyPvalue = methPercentageDiff(methdata.getValuesInRange(genedata[ensembl].chr,genedata[ensembl].start,genedata[ensembl].end))
            
            promotorClassification = classifyCpGContent(genome.getSequence(genedata[ensembl].chr,promotorstart,promotorstop).upper())
            geneBodyClassification = classifyCpGContent(genome.getSequence(genedata[ensembl].chr,genedata[ensembl].start,genedata[ensembl].end).upper())
            
            outputcsv.writerow([ensembl,genedata[ensembl].name,"1" if incpg else "0", affyavglogFC, affyavgfc, affyavgexpr, affyminp, promotorMethDiff, promotorPvalue, geneBodyMethDiff,genebodyPvalue,promotorClassification,geneBodyClassification])
        except sequence.genome.UnknownChromosomeException:
            continue