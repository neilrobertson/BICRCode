'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
from datastructures.AggregateTree import AggregateTree
import scipy.stats
from bed.treatment import SimpleBed
import csv
import sequence.genome
from sequence.genome import Genome


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["cpg-meth=","regions=","outputfile=","buildGenome=","threshold="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    methdata = None
    regionsfile= None
    outputfile = None
    genomeBuild = None
    threshold = None
    
    for o, a in opts:
        if o=="--cpg-meth":
            methdatafile = a
        elif o=="--regions":
            regionsfile = a
        elif o=="--outputfile":
            outputfile = a
        elif o=="--buildGenome":    
            genomeBuild = a
        #Integer as a threshold whether to print reads that have significantly low reads on right or left at this level   
        elif o=="--threshold":    
            threshold = a
            
    assert methdatafile != None
    assert regionsfile != None
    assert outputfile != None
    assert genomeBuild != None
    if threshold == None: 
        threshold = 10
    if isinstance(threshold, int) is False:
        threshold = 10  
    
    if genomeBuild not in ("hg19", "hg18", "mm9"): 
        genomeBuild = "hg19"
        print "Genome build type unacceptable. Defaulting to genome hg19..."
    
    regions = SimpleBed(regionsfile)
    methdata = AggregateTree(methdatafile)
    
    genome = Genome(genomeBuild)
    print "Completed loaded %s genome..." % (genomeBuild)
    
    def methChiSquared(pooledMeth):
        chi2, p = None, None
        try:
            chi2,p,dof,expected = scipy.stats.chi2_contingency(pooledMeth) #@UnusedVariable
        except ValueError:
            p = 1.0
        return chi2, p

    def methTotals(values):
        lmethTotal = 0
        lunmethTotal = 0
        rmethTotal = 0
        runmethTotal = 0
        
        for (lmeth,lunmeth,rmeth,runmeth) in values:
            lmethTotal += lmeth
            lunmethTotal += lunmeth
            rmethTotal += rmeth
            runmethTotal += runmeth
        
        return (lmethTotal,lunmethTotal,rmethTotal,runmethTotal)
    
    def methGetP(values):
        (lmeth,lunmeth,rmeth,runmeth) = methTotals(values)
        oddsratio, pvalue = scipy.stats.fisher_exact([[lmeth, lunmeth], [rmeth, runmeth]]) #@UnusedVariable
        #p = fisher_exact.pvalue(lmeth,lunmeth,rmeth,runmeth)
        return pvalue
        
        
    def methPercentageDiff(values, threshold):
        (lmeth,lunmeth,rmeth,runmeth) = methTotals(values)
        
        # arbitrary cutoff for regions with very few reads
        if lmeth + lunmeth <= threshold or rmeth + runmeth <= threshold:
            return None
        
        lmeth = float(lmeth)
        lunmeth = float(lunmeth)

        rmeth = float(rmeth)
        runmeth = float(runmeth)
        
        lpercentage = lmeth / (lmeth+lunmeth)
        rpercentage = rmeth / (rmeth+runmeth)
        return rpercentage - lpercentage
  
  
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
        Gs, Cs, cpgs, effectiveSeqLen, gcContent, cpgContent, cpgRatio = sequenceGCandCpGcontent(sequence) #@UnusedVariable
        
        hcp = False
        icp = False
        if gcContent >= 0.55 and cpgRatio >= 0.75:
            hcp = True
        elif cpgRatio >= 0.48:
            icp = True
            
        return "high" if hcp else ("intermediate" if icp else "low"), cpgs
  
  
    outputcsv = csv.writer(open(outputfile,"w"),delimiter="\t")
    outputcsv.writerow(["chr","start","stop","lmeth","lunmeth","rmeth","runmeth","methDiff","pValue","chiSquared"]) #,"classification","CpGCount"])
    
    for chrm,start,stop in regions:
        # look up meth value
        #print chrm, start, stop
        try:
            values = methdata.getValuesInRange(chrm,start,stop)
            
            (lmeth,lunmeth,rmeth,runmeth) = methTotals(values)
            
            methDiff = methPercentageDiff(values, threshold)
            if methDiff is not None:
                chi2,chi_P = methChiSquared(values)
                pValue = methGetP(values)
                seq = genome.getSequence(chrm,start,stop).upper()
                classification, cpgs = classifyCpGContent(seq)
                outputcsv.writerow([chrm,start,stop,lmeth,lunmeth,rmeth,runmeth,methDiff,pValue,chi_P,classification,cpgs])
        except sequence.genome.UnknownChromosomeException:
            continue