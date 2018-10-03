'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv
from itertools import izip_longest
from datastructures.memoized import memoized

# sudo easy_install mpmath
import mpmath
#mpmath.mp.dps = 20 # decimal digits of precision

@memoized
def mpmathcdf(g,df,dps=10):
    mpmath.mp.dps = dps
    
    x,k = mpmath.mpf(g), mpmath.mpf(df)
    cdf = mpmath.gammainc(k/2, 0, x/2, regularized=True)
    
    # floating point precision insufficient, use more precision
    if cdf == 1.0:
        if dps > 4000:
            return cdf # give up after a certain point
        else:
            cdf = mpmathcdf(g,df,dps*2)
    return cdf

def GtoP(g,df):
    assert g >= 0, g
    return float(1-mpmathcdf(g,df))

import math

@memoized
def flnf(f):
    return f*math.log(f) if f>0.5 else 0

@memoized
def gtest(a,b,c,d):
    row1 = a+b
    row2 = c+d
    col1 = a+c
    col2 = b+d
    
    total = flnf(a+b+c+d)
    celltotals = flnf(a)+flnf(b)+flnf(c)+flnf(d)
    rowtotals = flnf(row1)+flnf(row2)
    coltotals = flnf(col1)+flnf(col2)
    
    # abs is necessary to account for float precision errors when doing the subtraction
    # in some cases where rowtotals and coltotals are the same the float subtraction
    # 
    return abs(2*((celltotals + total)-(rowtotals+coltotals)))


# total G: at least some of the data is different
# 1 degree of freedom per replicate (i.e. 3 in total), add G's from each seperate test
def totalG(data):
    dof = 0
    totalg = 0
    for lmeth,lunmeth,rmeth,runmeth in grouper(4,data): # extract in quads (2 pairs)
            totalg += gtest(lmeth,lunmeth,rmeth,runmeth)
            dof+=1 
    assert totalg >= 0
    return totalg, dof, GtoP(totalg,dof)

# pooled G: the data as a whole is different
# 1 degree of freedom, pool results from all replicates and do 1 2x2 table
def pooledG(data):
    lmethtotal = 0
    lunmethtotal = 0
    rmethtotal = 0
    runmethtotal = 0
    
    for lmeth,lunmeth,rmeth,runmeth in grouper(4,data): # extract in quads (2 pairs)
        lmethtotal += lmeth
        lunmethtotal += lunmeth
        rmethtotal += rmeth
        runmethtotal += runmeth
    
    dof = 1
    
    pooledg = gtest(lmethtotal,lunmethtotal,rmethtotal,runmethtotal)
    
    assert pooledg >= 0
    
    return pooledg, dof, GtoP(pooledg,dof)

# there is evidence that they are significantly different from each other
# heterogeneity G: total G - pooled G, df = 3-1 = 2
def heterogeneityG(totalg,totaldof,pooledg,pooleddof):
    heteroG = max(totalg-pooledg,0)
    heterodof = totaldof-pooleddof
    return heteroG,heterodof,GtoP(heteroG,totaldof-pooleddof)



def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def unpackLine(chrm,coord,*data):
    return chrm.strip(), int(coord), [int(d) for d in data]

def extractLine(line):
    chrm,coord,data = unpackLine(*line)
    if not chrm.startswith("chr"):
        chrm = "chr"+chrm
    return (chrm,coord,data)


    
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
    
    
    minReadsInBothSamples = 3
    minReadsInOneSample = 10
    depthLowerRatio = 0.8
    depthUpperRatio = 1.2
    
    minWindowSize = 1000
    maxWindowSize = 5000
    maxSkipSize = 2000
    alwaysExtendDistance = 10
    
    currentWindowPos = None
    currentChr = None
    
    cs = 0
    dataCount = None
    
    def totalMethCounts(data):
        lmethtotal = 0
        lunmethtotal = 0
        
        rmethtotal = 0
        runmethtotal = 0
        
        for lmeth,lunmeth,rmeth,runmeth in grouper(4,data): # extract in quads (2 pairs)
            lmethtotal += lmeth
            lunmethtotal += lunmeth
            rmethtotal += rmeth
            runmethtotal += runmeth
        
        return lmethtotal,lunmethtotal,rmethtotal,runmethtotal
    
    def writeWindow(chrm,pos,epos):        
        global dataCount
        global cs
        global windows
        
        lmethCount,lunmethCount,rmethCount,runmethCount = totalMethCounts(dataCount)
        
        if (lmethCount + lunmethCount) == 0 or (rmethCount + runmethCount) == 0:
            print "Unprintable window:",chrm, pos, epos, lmethCount, lunmethCount, rmethCount, runmethCount
            return
        
        lPercentage = (float(lmethCount)/float(lmethCount+lunmethCount)) * 100.0
        rPercentage = (float(rmethCount)/float(rmethCount+runmethCount)) * 100.0
        
        total,totaldof,totalp = totalG(dataCount)
        pooled,pooleddof,pooledp = pooledG(dataCount)
        hetero,heterodof,heterop = heterogeneityG(total,totaldof,pooled,pooleddof)
        
        windows.writerow([chrm,pos,epos,'1' if rPercentage > lPercentage else '-1',totalp,pooledp,heterop,cs])
        
        dataCount = None
        cs = 0

    lastCoordCounted = None
    # the actual work here
    for line in methfile:
        
        chrm,coord,data = extractLine(line)

        if currentChr == None:
            currentChr = chrm
            currentWindowPos = coord
            lastCoordCounted = coord
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
            elif coord - currentWindowPos > maxWindowSize:
                writeWindow(currentChr,currentWindowPos,lastCoordCounted)
                currentChr = chrm
                currentWindowPos = coord
            
            # otherwise only consider writing a window if it's big enough
            elif lastCoordCounted - currentWindowPos > minWindowSize and cs > 6:
                # next coord must not be right next to the window or we should just continue it
                if coord-lastCoordCounted > alwaysExtendDistance:
                    # and has enough depth in both samples
                    lmeth,lunmeth,rmeth,runmeth = totalMethCounts(dataCount)
                    
                    if (lmeth + lunmeth >= minReadsInBothSamples) and (rmeth + runmeth >= minReadsInBothSamples):
                        # and enough depth in at least one sample
                        if (lmeth + lunmeth > minReadsInOneSample or rmeth + runmeth > minReadsInOneSample):
                            # and depth is similar between the two samples
                            depthRatio = float(lmeth+lunmeth)/float(rmeth+runmeth)
                            if depthRatio >= depthLowerRatio and depthRatio <= depthUpperRatio:
                                writeWindow(currentChr,currentWindowPos,lastCoordCounted)
                                currentChr = chrm
                                currentWindowPos = coord
        
        ####
        # add data to existing data count
        ####
        if dataCount == None:
            dataCount = data
        else:
            for index,d in enumerate(data):
                dataCount[index] += d
        lastCoordCounted = coord
        
        # inc number of c's counted
        cs += 1
        
    # after we write a window there is always at least one unmeasured coord after it so we have one more window to write
    writeWindow(currentChr,currentWindowPos,lastCoordCounted)