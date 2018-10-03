'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv

# https://github.com/brentp/fishers_exact_test 
from fisher import cfisher as fisher_exact

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
    
    p = fisher_exact.pvalue(lmeth,lunmeth,rmeth,runmeth).two_tail
    #return -10.0*math.log10(p)
    return p

    
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["combined="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    
    for o, a in opts:
        if o=="--combined":
            infile = a
            print "combined", a
    
    assert infile != None
    
    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    
    # ---    
    # windows
    
    windows = csv.writer(open(infile+".windows.bedGraph","w"),delimiter='\t')
    
    
    minReadsInBothSamplesForFisher = 1
    minReadsInOneSampleForFisher = 3
    fisherDepthLowerRatio = 0.8
    fisherDepthUpperRatio = 1.2
    
    minWindowSize = 200
    maxSkipSize = 2000
    
    currentWindowPos = None
    currentChr = None
    
    lmethCount = 0
    lunmethCount = 0
    rmethCount = 0
    runmethCount = 0
    cs = 0
    
    def writeWindow(chrm,pos,epos):
        
        global lmethCount
        global lunmethCount
        global rmethCount
        global runmethCount
        global cs
        
        if (lmethCount + lunmethCount) == 0 or (rmethCount + runmethCount) == 0:
            #print "Unprintable window:",chrm, pos, epos, lmethCount, lunmethCount, rmethCount, runmethCount
            return
        
        lPercentage = (float(lmethCount)/float(lmethCount+lunmethCount)) * 100.0
        rPercentage = (float(rmethCount)/float(rmethCount+runmethCount)) * 100.0
        
        windows.writerow([chrm,pos,epos,'1' if rPercentage > lPercentage else '-1',fisherExact(lmethCount,lunmethCount,rmethCount,runmethCount),cs])
        
        lmethCount = 0
        lunmethCount = 0
        rmethCount = 0
        runmethCount = 0
        cs = 0

    lastCoordCounted = None
    # the actual work here
    for line in methfile:
        #print line
        (chrm,coord,lmeth,lunmeth,rmeth,runmeth) = parseMethLine(line)

        if currentChr == None:
            currentChr = chrm
            currentWindowPos = coord
        else:
            # always write a window if we change chrm
            if chrm != currentChr:
                writeWindow(currentChr,currentWindowPos,lastCoordCounted)
                currentChr = chrm
                currentWindowPos = coord
            # always write a window if there is too much distance to the next c to be counted (ie a big gap)
            elif coord - lastCoordCounted > maxSkipSize:
                writeWindow(currentChr,currentWindowPos,lastCoordCounted)
                currentChr = chrm
                currentWindowPos = coord
                
            # always write a window if the current window would get too big
            #elif coord - currentWindowPos > maxWindowSize:
            #    writeWindow(currentChr,currentWindowPos,lastCoordCounted)
            #    currentChr = chrm
            #    currentWindowPos = coord
            
            # otherwise only consider writing a window if it's big enough
            elif currentWindowPos+minWindowSize < coord:
                # and has enough depth in both samples
                if (lmeth + lunmeth >= minReadsInBothSamplesForFisher) and (rmeth + runmeth >= minReadsInBothSamplesForFisher):
                    # and enough depth in at least one sample
                    if (lmeth + lunmeth > minReadsInOneSampleForFisher or rmeth + runmeth > minReadsInOneSampleForFisher):
                        # and depth is similar between the two samples
                        depthRatio = float(lmeth+lunmeth)/float(rmeth+runmeth)
                        if depthRatio >= fisherDepthLowerRatio and depthRatio <= fisherDepthUpperRatio:
                                writeWindow(currentChr,currentWindowPos,lastCoordCounted) # half open coords, ie. coord isn't included in the window
                                currentChr = chrm
                                currentWindowPos = coord
        
        lastCoordCounted = coord
                                  
        lmethCount += lmeth
        lunmethCount += lunmeth
        rmethCount += rmeth
        runmethCount += runmeth
        cs += 1
        
    # after we write a window there is always at least one unmeasured coord after it so we have one more window to write
    writeWindow(currentChr,currentWindowPos,lastCoordCounted)