'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
from datastructures.AggregateTree import AggregateTree
from csvfile.indexedcsv import IndexedCSV
from genemapping.Ensembl import EnsemblGenes
import csv
import os
from bed.treatment import ExtendedBed
import sequence.genome
from sequence.genome import Genome
from csvfile.genelist import GeneList
from itertools import izip


# https://github.com/brentp/fishers_exact_test 
from fisher import cfisher as fisher_exact

def fisherExact(lmeth,lunmeth,rmeth,runmeth):
    
    p = fisher_exact.pvalue(lmeth,lunmeth,rmeth,runmeth).two_tail
    #return -10.0*math.log10(p)
    return p

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["cpg-meth-file=",
                                                      "gene-expression-file=",
                                                      "fccol=",
                                                      "exprcol=",
                                                      "outputfile=",
                                                      "promotorsize=",
                                                      "reverseOrder",
                                                      "genelist=",
                                                      "assembly="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    methdatafile = None
    
    exprfile= None
    ensemblidcol = "ensemblid"
    sigcol = "q_value"
    fccol = "log2(fold_change)"
    exprcols = ["value_1","value_2"]
        
    outputfile = None
    
    upstreamPromotor = 5000
    downstreamPromotor = 1000
    
    tsssize = 500
    
    reverse = False
    
    genelists = []
    
    assembly = "hg18"
    
    for o, a in opts:
        if o=="--cpg-meth-file":
            methdatafile = a
        elif o=="--gene-expression-file":
            exprfile = a
        elif o=="--fccol":
            fccol = a
        elif o=="--exprcols":
            exprcols = a.split(",")
        elif o=="--outputfile":
            outputfile = a
        elif o=="--promotorsize":
            upstreamPromotor = int(a)
            downstreamPromotor = int(a)
        elif o=="--reverseOrder":
            reverse = True
        elif o=="--genelist":
            genelists.append(GeneList(a))
        elif o=="--assembly":
            assembly = a
    
    assert methdatafile != None
    
    # if we have an expression file we need fc and expression columns
    if exprfile != None:
        assert fccol != None
        assert exprcols != None
    else:
        exprcols = []
    
    assert outputfile != None
    
    genedata = EnsemblGenes(assembly=assembly)
    
    genome = Genome(genomeBuild = assembly)
    
    if assembly == "hg18":        
        cpgIslands = ExtendedBed(os.path.expanduser("~/mount/publicdata/hg18/cpgislands/cpgislands.bed"))
        lINEs = ExtendedBed(os.path.expanduser("~/mount/publicdata/hg18/repeats/LINEs-0.bed"))
        sINEs = ExtendedBed(os.path.expanduser("~/mount/publicdata/hg18/repeats/SINEs-0.bed"))
    elif assembly == "hg19":
        cpgIslands = ExtendedBed(os.path.expanduser("~/mount/publicdata/hg19/CpGIslands/cpgislands.bed"))
        lINEs = ExtendedBed(os.path.expanduser("~/mount/publicdata/hg19/Repeats/UCSC_HG19_LINEs.bed"),defaultkeys=["chrom","chromStart","chromEnd","name","strand"])
        sINEs = ExtendedBed(os.path.expanduser("~/mount/publicdata/hg19/Repeats/UCSC_HG19_SINEs.bed"),defaultkeys=["chrom","chromStart","chromEnd","name","strand"])
    else:
        assert False, "Unknown genome build for cpgislands / lines / sines"
    
    
    methdata = AggregateTree(methdatafile, reverse = reverse)
    
    def methPercentageDiff(values):
        
        assert values != None
        
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
        
        rmeth = float(rmethTotal)
        runmeth = float(runmethTotal)
        
        # arbitrary cutoff for regions with very few reads
        if lmeth + lunmeth <= 10.0 or rmeth + runmeth <= 10.0:
            return None,None,None,None,None
        
        lpercentage = lmeth / (lmeth+lunmeth)
        rpercentage = rmeth / (rmeth+runmeth)
        
        pvalue = fisherExact(lmethTotal,lunmethTotal,rmethTotal,runmethTotal)
        
        # add a small psuedocount to both numerator and denominator to avoid div by 0 for relative meth
        return lpercentage,rpercentage,((rpercentage - lpercentage)+0.01)/(lpercentage+0.01), rpercentage - lpercentage, pvalue

  
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
            
        return "high" if hcp else ("intermediate" if icp else "low"), cpgRatio
    
    #Method to determine the proportion of a feature type (chr:start-stop) which is occupied by another feature type (extendedBed)
    def calssifyIntervalContent(chr, start, end, extendedBed):
        intervals = extendedBed.getIntervalsInRange(chr, start, end)
        density = 0
        if (len(intervals) > 0): 
            for interval in intervals:
                if (interval.start < start):
                    density += interval.end - start
                if (interval.end > end):
                    density += end - interval.start
                if (interval.start > start and interval.end < end):
                    density += interval.end - interval.start
        return float(density)/float(end - start)

  

  
    outputcsv = csv.writer(open(outputfile,"w"),delimiter="\t")
    
    header = ["ensemblid","genename","cpgisland"]
    
    if exprfile != None:
        header.extend(["logfc","sig"])
        header.extend(exprcols)
        
    header.extend(["promotorMethDiff","relativePromotorMethDiff","promotorPvalue",
                   "genebodyMethDiff","relativeGenebodyMethDiff","genebodyPvalue",
                   "tssMethDiff","relativeTssMethDiff","tssPvalue",
                   "promotorCpGClassification","genebodyCpGClassification","tssCpGClassification"])
    header.extend(["lpromotorMeth","rpromotorMeth",
                   "lgenebodyMeth","rgenebodyMeth",
                   "ltssMeth","rtssMeth"])
    header.extend(["promotorCpGRatio","genebodyCpGRatio","tssCpGRatio"])
    header.extend(["nearestIsland", "islandStart", "islandEnd", "islandLength", "islandlMeth","islandrMeth"])
    header.extend(["shoresl", "shoresr", "shelvesl", "shelvesr"])
    header.extend(["lineContentGene", "sineContentGene", "lineContentPromoter", "sineContentPromoter"])
    
    for genelist in genelists:
        header.append(genelist.getFriendlyName())
    
    outputcsv.writerow(header)
    
    # if expression file then use those genes for iterator
    # else use all ensembl genes
    if exprfile != None:
        exprCSV = IndexedCSV(exprfile,keyPos=1)
        iter = exprCSV
    else:
        iter = genedata 
    
    for testid in iter:
        # look up meth value
        
        try:
            
            if exprfile != None:
                ensembl = exprCSV[testid][ensemblidcol]
                
                if ensembl not in genedata:
                    # RNA-Seq is not a single unique gene, skip
                    continue
                
                logFC = exprCSV[testid][fccol]
                sig = float(exprCSV[testid][sigcol])<=0.05
            else:
                ensembl = testid
            

            
            promotorstart, promotorstop = genedata[ensembl].getPromotorRegion(upstreamPadding = upstreamPromotor, downstreamPadding = downstreamPromotor)
            
            incpg = len(cpgIslands.getValuesInRange(genedata[ensembl].chr, promotorstart, promotorstop))!=0    
            
            #Nearest Island
            islandlmeth,islandrmeth,islandStart,islandEnd,islandLength   = "","","","",""
            nearestCpGIsland, nearestCpGIslandObject = None,None
            nearbyIslands = cpgIslands.getIntervalsInRange(genedata[ensembl].chr, genedata[ensembl].start-5000, genedata[ensembl].start+5000)
            
            if (len(nearbyIslands) > 0): 
                for island in nearbyIslands:
                    
                    distance = 0
                    if (island.start <= genedata[ensembl].start and island.end >= genedata[ensembl].start):
                        distance = 0
                    elif (island.start >= genedata[ensembl].start):
                        distance = island.start-genedata[ensembl].start    
                    elif (island.end <= genedata[ensembl].start):
                        distance = genedata[ensembl].start - island.end
                        
                    if (distance < nearestCpGIsland or nearestCpGIsland == None):
                        nearestCpGIsland,nearestCpGIslandObject = distance,island
            
            if (nearestCpGIslandObject != None):
                islandlmeth, islandrmeth,_,_,_ = methPercentageDiff(methdata.getValuesInRange(genedata[ensembl].chr,nearestCpGIslandObject.start,nearestCpGIslandObject.end))
                islandStart,islandEnd,islandLength = nearestCpGIslandObject.start,nearestCpGIslandObject.end, nearestCpGIslandObject.end-nearestCpGIslandObject.start
                
            #SINE & LINE density - make into a helper method
            lineContentGene = calssifyIntervalContent(genedata[ensembl].chr, genedata[ensembl].start, genedata[ensembl].end, lINEs)
            sineContentGene = calssifyIntervalContent(genedata[ensembl].chr, genedata[ensembl].start, genedata[ensembl].end, sINEs)
            lineContentPromotor = calssifyIntervalContent(genedata[ensembl].chr,promotorstart,promotorstop, lINEs)
            sineContentPromotor = calssifyIntervalContent(genedata[ensembl].chr,promotorstart,promotorstop, sINEs)
              
            #Shelves and Shores - Holger Heyn,Manel Esteller, DNA methylation profiling in the clinic:applications and challenges, Nature Genetics Reviews, 2012
            shoresl, shoresr, shelvesl, shelvesr = None,None,None,None
            if (nearestCpGIslandObject != None):
                try:
                    shoresl, shoresr,_,_,_ =  [(i+j)/2 for i,j in izip(methPercentageDiff(methdata.getValuesInRange(genedata[ensembl].chr,nearestCpGIslandObject.start-2000,nearestCpGIslandObject.start)),methPercentageDiff(methdata.getValuesInRange(genedata[ensembl].chr,nearestCpGIslandObject.end,nearestCpGIslandObject.end+2000)))]
                    shelvesl, shelvesr,_,_,_ = [(i+j)/2 for i,j in izip(methPercentageDiff(methdata.getValuesInRange(genedata[ensembl].chr,nearestCpGIslandObject.start-4000,nearestCpGIslandObject.start-2000)),methPercentageDiff(methdata.getValuesInRange(genedata[ensembl].chr,nearestCpGIslandObject.end+2000,nearestCpGIslandObject.end+4000)))]
                except TypeError:
                    print 'Shelve/Shore out of bounds' #No value returned?
                
            
            lpromotorMeth,rpromotorMeth,relativePromotorMethDiff,promotorMethDiff,promotorPvalue = methPercentageDiff(methdata.getValuesInRange(genedata[ensembl].chr,promotorstart,promotorstop))
            lgenebodyMeth,rgenebodyMeth,relativeGeneBodyMethDiff,geneBodyMethDiff,geneBodyPvalue = methPercentageDiff(methdata.getValuesInRange(genedata[ensembl].chr,genedata[ensembl].start,genedata[ensembl].end))
            ltssMeth,     rtssMeth,     relativeTssMethDiff,     tssMethDiff,     tssPvalue      = methPercentageDiff(methdata.getValuesInRange(genedata[ensembl].chr,genedata[ensembl].start-tsssize/2,genedata[ensembl].start+tsssize/2))
            
            promotorClassification, promotorCpGRatio = classifyCpGContent(genome.getSequence(genedata[ensembl].chr,promotorstart,promotorstop).upper())
            geneBodyClassification, genebodyCpGRatio = classifyCpGContent(genome.getSequence(genedata[ensembl].chr,genedata[ensembl].start,genedata[ensembl].end).upper())
            tssClassification,      tssCpGRatio      = classifyCpGContent(genome.getSequence(genedata[ensembl].chr,genedata[ensembl].start-tsssize/2,genedata[ensembl].start+tsssize/2).upper())
            
            
            row =  [ensembl,genedata[ensembl].name,"1" if incpg else "0"]
            
            if exprfile != None:
                row.extend([logFC, "1" if sig else "0"])
            
            exprs = []
            for exprcol in exprcols:
                exprs.append(exprCSV[testid][exprcol])
            row.extend(exprs)
            
            row.extend([promotorMethDiff,relativePromotorMethDiff,promotorPvalue,
                        geneBodyMethDiff,relativeGeneBodyMethDiff,geneBodyPvalue,
                        tssMethDiff,relativeTssMethDiff,tssPvalue,
                        promotorClassification,geneBodyClassification,tssClassification])
            
            row.extend([lpromotorMeth,rpromotorMeth,
                        lgenebodyMeth,rgenebodyMeth,
                        ltssMeth,rtssMeth])
            
            row.extend([promotorCpGRatio,genebodyCpGRatio,tssCpGRatio]) 
            row.extend([nearestCpGIsland,islandStart,islandEnd,islandLength, islandlmeth, islandrmeth])
            row.extend([shoresl,shoresr,shelvesl,shelvesr])
            row.extend([lineContentGene,sineContentGene,lineContentPromotor,sineContentPromotor])
            
            for genelist in genelists:
                if ensembl in genelist.seengenes:
                    row.append("1")
                else:
                    row.append("0")
            
            outputcsv.writerow(row)
        except sequence.genome.UnknownChromosomeException:
            continue