'''
Created on 25 Jul 2014

@author: neilrobertson
'''

import getopt
import sys
import scipy.stats
from datastructures.PooledAggregateTree import PooledAggregateTree
from bed.treatment import SimpleBed
import csv
import sequence.genome
from sequence.genome import Genome


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["cpg-meth=","regions=","outputfile=","buildGenome=","replicateAnnotation=","replicateTypes="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    methdata = None
    regionsfile= None
    outputfile = None
    genomeBuild = None
    replicateAnnotation = None
    replicateTypes = None
    
    for o, a in opts:
        if o=="--cpg-meth":
            methdatafile = a
        elif o=="--regions":
            regionsfile = a
        elif o=="--outputfile":
            outputfile = a
        elif o=="--buildGenome":    
            genomeBuild = a
        elif o=="--replicateAnnotation":
            replicateAnnotation = a
        elif o=="--replicateTypes":
            replicateTypes = a
            
    assert methdatafile != None
    assert regionsfile != None
    assert outputfile != None
    assert genomeBuild != None
    assert replicateAnnotation != None
    assert replicateTypes != None
    
    
    if len(replicateAnnotation)%len(replicateTypes) == 0:
        print "Unbalanced length of replicant annotators..."    
    if (len(replicateAnnotation)/2) != len(replicateTypes):
        print "Unbalanced counts of replicant annotators to replicant types..."    
    if genomeBuild not in ("hg19", "hg18", "mm9"): 
        genomeBuild = "hg19"
        print "Genome build type unacceptable. Defaulting to genome hg19..."
    
    print "Loading regions file into memory..."
    regions = SimpleBed(regionsfile)
    print "Loading methylation data into memory..."
    methdata = PooledAggregateTree(methdatafile)
    print "Loading genome into structure..."
    genome = Genome(genomeBuild)
    
    
    def methGetP(rep1, rep2):
        ''' 
        Calculate fishers exact test (p-value) from the methylated and unmethylated 
        scores between two sets of replicates
        '''
        oddsRatio, p = scipy.stats.fisher_exact([rep1, rep2]) #@UnusedVariable
        return p
        

    def methPercentageDiff(rep1, rep2):
        ''' 
        Calculate the percentage difference between 
        two replicates percentage methylation scores
        '''
        if rep1[0] + rep1[1] <= 10 or rep2[0] + rep2[1] <= 10:
            return " "
        
        lmeth = float(rep1[0])
        lunmeth = float(rep1[1])
        rmeth = float(rep2[0])
        runmeth = float(rep2[1])
        
        lpercentage = lmeth / (lmeth+lunmeth)
        rpercentage = rmeth / (rmeth+runmeth)
        return rpercentage - lpercentage
    
    
    
    def methChiSquared(pooledMeth):
        chi2, p = None, None
        try:
            chi2,p,dof,expected = scipy.stats.chi2_contingency(pooledMeth) 
        except ValueError:
            p = 1.0
        return chi2, p
    
        
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
        return "high" if hcp else ("intermediate" if icp else "low"), cpgs
  
  
  
    def buildReplicantsDict(replicateAnnotation, replicateTypes):
        typeDict = {}
        finalReplicateAnnotators = []
        for type in replicateAnnotation.split(","):
            count = 1
            try:
                count = typeDict[type]
                count = count + 1
                typeDict[type] = count
            except:
                typeDict[type] = count
            finalReplicateAnnotators.append("%s%s" % (type, str(count)))
        return finalReplicateAnnotators
    
    
    def pooledMethTotals(values, resultsDict, finalReplicateAnnotators):
        
        for uniqueRead in values:
            counter = 0
            repCounter = 0
            repMethylation = []
            for i in uniqueRead:
                if counter%2 == 0:
                    repMethylation.append(i)
                else:
                    repMethylation.append(i)
                    meth = int(repMethylation[0])
                    unmeth = int(repMethylation[1])
                    replicate = finalReplicateAnnotators[repCounter]
                    currentTotalMeth = resultsDict[replicate][0] + meth
                    currentTotalUnmeth = resultsDict[replicate][1] + unmeth
                    resultsDict[replicate] = [currentTotalMeth, currentTotalUnmeth]
                    repCounter += 1
                    repMethylation = []
                    if repCounter == len(finalReplicateAnnotators):
                        repCounter = 0
                counter += 1   
        return resultsDict   
     
     
     
    def testing(values):
        test = [0,0,0,0,0,0,0,0]
        for uniqueRead in values:
            count = 0
            for i in uniqueRead:
                current = test[count]
                current += int(i)
                test[count] = current
                count += 1
        return test
    
    
    finalReplicateAnnotators = buildReplicantsDict(replicateAnnotation, replicateTypes)
    
    #Create headers for the combined CSV file             
    outputcsv = csv.writer(open(outputfile,"w"),delimiter="\t")
    headers = ["chr","start","stop"]
    for a in finalReplicateAnnotators:
        headers.append("%s.meth" % (a))
        headers.append("%s.unmeth" % (a))
    for c in replicateTypes.split(","):
        headers.append("%s.pooledMeth" % (c))
        headers.append("%s.pooledUnmeth" % (c)) 
    for s in replicateTypes.split(","):
        headers.append("%s.chi_Squared" % (s))   
    for s in replicateTypes.split(","):
        for i in replicateTypes.split(","):
            if s != i:
                headers.append("%s.%s_p_value" % (s,i))
                headers.append("%s.%s_methDiff" % (s,i))
    for b in finalReplicateAnnotators:
        for c in finalReplicateAnnotators:
            if b != c:
                headers.append("%s.%s_p_value" % (b,c))    
                headers.append("%s.%s_methDiff"% (b,c))  
    headers.append("cpgs")
    outputcsv.writerow(headers)
    
    currentChrm = None
    
    for chrm,start,stop in regions:
        if chrm != currentChrm:
            currentChrm = chrm
            print "Working on %s" % (currentChrm)
            
        resultsDict = {}
        for i in finalReplicateAnnotators:
            resultsDict[i] = [0,0] #meth/unmeth
            
        try:
            values = methdata.getValuesInRange(chrm,start,stop)
            
            if len(values) > 0:
                output = [chrm, start, stop]
                #test = testing(values)
                #print test
                resultsDict = pooledMethTotals(values, resultsDict, finalReplicateAnnotators)
                for s in finalReplicateAnnotators:
                    output.append(resultsDict[s][0])
                    output.append(resultsDict[s][1])
                    
                pooledTotals = {}
                for replicateType in replicateTypes.split(","):
                    for t in finalReplicateAnnotators:
                        if t.find(replicateType) != -1:
                            try:
                                counts = pooledTotals[replicateType]
                                counts = [counts[0] + resultsDict[t][0], counts[1] + resultsDict[t][1]]
                                pooledTotals[replicateType] = counts
                            except:
                                pooledTotals[replicateType] = [resultsDict[t][0], resultsDict[t][1]]   
                                finalReplicateAnnotators
                for s in replicateTypes.split(","):
                    output.append(pooledTotals[s][0])
                    output.append(pooledTotals[s][1])
                    
                chi2Tots = {}    
                for s in replicateTypes.split(","):
                    for i in finalReplicateAnnotators:
                        if i.find(s) != -1:
                            try:
                                chi2Tots[s][0].append(resultsDict[i][0])
                                chi2Tots[s][1].append(resultsDict[i][1])
                            except:
                                chi2Tots[s] = [[resultsDict[i][0]], [resultsDict[i][1]]]
                  
                for s in replicateTypes.split(","):
                    pooledMeths = chi2Tots[s]
                    chi2, p = methChiSquared(pooledMeths)
                    output.append(p)

                for s in replicateTypes.split(","):
                    for i in replicateTypes.split(","):
                        if s != i:
                            output.append(str(methGetP(pooledTotals[s], pooledTotals[i])))  
                            output.append(str(methPercentageDiff(pooledTotals[s], pooledTotals[i]))) 
                             
                for x in finalReplicateAnnotators:
                    firstRep = resultsDict[x]
                    for y in finalReplicateAnnotators:
                        if y != x:
                            secondRep = resultsDict[y]
                            output.append(str(methGetP(firstRep, secondRep)))  
                            output.append(str(methPercentageDiff(firstRep, secondRep)))  

                seq = genome.getSequence(chrm,start,stop).upper()
                classification, cpgs = classifyCpGContent(seq)
                output.append(cpgs)
                outputcsv.writerow(output)
        except sequence.genome.UnknownChromosomeException:
            continue