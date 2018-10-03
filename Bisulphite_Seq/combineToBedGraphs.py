
'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv
import math
import scipy.stats as sci_stat
# https://github.com/brentp/fishers_exact_test 
#from fisher import cfisher as fisher_exact

def parseMethLine(line):
    (chrm,coord,lmeth,lunmeth,rmeth,runmeth) = line
    chrm = chrm.strip()
    if not chrm.startswith("chr"):
        chrm = "chr"+chrm
    # coord should be in base 0 already, just need it as an int
    coord = int(coord)
    lmeth=int(lmeth)
    lunmeth=int(lunmeth)
    rmeth=int(rmeth)
    runmeth=int(runmeth)
    return (chrm,coord,lmeth,lunmeth,rmeth,runmeth)


def fisherExact(lmeth,lunmeth,rmeth,runmeth):
    
    odds, p = sci_stat.fisher_exact([[lmeth,lunmeth],[rmeth,runmeth]], alternative="two-sided") 
    #return -10.0*math.log10(p)
    return p

    
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["combined=","leftsuffix=","rightsuffix="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    leftsuffix = None
    rightsuffix = None 
    
    for o, a in opts:
        if o=="--combined":
            infile = a
            print "combined", a
        elif o=="--leftsuffix":
            leftsuffix = a
            print "Left suffix", a
        elif o=="--rightsuffix":
            rightsuffix = a
            print "Right suffix", a

    
    assert infile != None
    assert leftsuffix != None
    assert rightsuffix != None
    
    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    
    # ---
    
    # for each
    # coverage
    # percentage
    # meth score
    
    # difference
    
    leftMethScore = csv.writer(open(infile+"-"+leftsuffix+".meth.bedGraph","w"),delimiter='\t')
    leftPercentage = csv.writer(open(infile+"-"+leftsuffix+".percent.bedGraph","w"),delimiter='\t')
    leftCoverage = csv.writer(open(infile+"-"+leftsuffix+".coverage.bedGraph","w"),delimiter='\t')

    rightMethScore = csv.writer(open(infile+"-"+rightsuffix+".meth.bedGraph","w"),delimiter='\t')
    rightPercentage = csv.writer(open(infile+"-"+rightsuffix+".percent.bedGraph","w"),delimiter='\t')
    rightCoverage = csv.writer(open(infile+"-"+rightsuffix+".coverage.bedGraph","w"),delimiter='\t')
    
    fisherP = csv.writer(open(infile+".fisher.bedGraph","w"),delimiter='\t')
    
    minReadsInAtLeastOneSampleForPercentage = 10
    minReadsInBothSamplesForPercentage = 3
    
    minReadsInBothSamplesForFisher = 1
    minReadsInOneSampleForFisher = 3
    fisherDepthLowerRatio = 0.8
    fisherDepthUpperRatio = 1.2

    
    def methScore(methylated,unmethylated):
        methScore = methylated - unmethylated
        # log the methScore (note you can't log a 0 or a negative number so do some fixing for that
        if methScore != 0:
            if methScore > 0:
                methScore = math.log(methScore)
            else:
                methScore = -1*math.log(abs(methScore))
        return methScore
    

    for line in methfile:
        (chrm,coord,lmeth,lunmeth,rmeth,runmeth) = parseMethLine(line)
        
        leftMethScore.writerow([chrm,coord,coord+1,methScore(lmeth,lunmeth)])
        rightMethScore.writerow([chrm,coord,coord+1,methScore(rmeth,runmeth)])        

        
        # only write the percentages if at least minReads number of reads
        if (lmeth + lunmeth >= minReadsInAtLeastOneSampleForPercentage or rmeth + runmeth >= minReadsInAtLeastOneSampleForPercentage) and (lmeth + lunmeth >= minReadsInBothSamplesForPercentage and rmeth + runmeth >= minReadsInBothSamplesForPercentage):
            
            lPercentage = (float(lmeth)/float(lmeth+lunmeth)) * 100.0
            rPercentage = (float(rmeth)/float(rmeth+runmeth)) * 100.0
            
            leftPercentage.writerow([chrm, coord, coord+1, lPercentage])
            rightPercentage.writerow([chrm, coord, coord+1, rPercentage])
        
        leftCoverage.writerow([chrm,coord,coord+1,lmeth+lunmeth])
        rightCoverage.writerow([chrm,coord,coord+1,rmeth+runmeth])
        
        if (lmeth + lunmeth >= minReadsInBothSamplesForFisher) and (rmeth + runmeth >= minReadsInBothSamplesForFisher):
            depthRatio = float(lmeth+lunmeth)/float(rmeth+runmeth)
            
            lPercentage = (float(lmeth)/float(lmeth+lunmeth)) * 100.0
            rPercentage = (float(rmeth)/float(rmeth+runmeth)) * 100.0
            
            if depthRatio >= fisherDepthLowerRatio and depthRatio <= fisherDepthUpperRatio:
                if (lmeth + lunmeth > minReadsInOneSampleForFisher or rmeth + runmeth > minReadsInOneSampleForFisher):
                    fisherP.writerow([chrm,coord,coord+1,'1' if rPercentage > lPercentage else '-1', fisherExact(lmeth,lunmeth,rmeth,runmeth)])