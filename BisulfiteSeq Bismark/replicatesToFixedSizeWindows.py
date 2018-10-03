'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv
from itertools import izip_longest
from itertools import izip
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



# extract the pooled meth counts from the sample
def pooledMethCounts(data,covs,compare):
    lmethtotal = 0
    lunmethtotal = 0
    rmethtotal = 0
    runmethtotal = 0
    
    if compare != None and len(compare)==2:
        # if we are given covariates to use
        covleft, covright = compare    
        for cov, (meth,unmeth) in izip(covs,grouper(2,data)): # extract in pairs of meth/unmeth
            if cov == covleft:
                lmethtotal += meth
                lunmethtotal += unmeth
            elif cov == covright:
                rmethtotal += meth
                runmethtotal += unmeth
    else:
        # if we have not been given them, assume that each adjacent sample (which is one pair) is one replicate
        for lmeth,lunmeth,rmeth,runmeth in grouper(4,data): # extract in quads (2 pairs) representing one replicate
            lmethtotal += lmeth
            lunmethtotal += lunmeth
            rmethtotal += rmeth
            runmethtotal += runmeth
        
    return lmethtotal,lunmethtotal,rmethtotal,runmethtotal


# extract individual replicate counts from the sample and return as a list of (lmeth,lunmeth,rmeth,runmeth)
# where l and r refer to the left and right cov in the compare array
# if compare array is empty assume that data is formatted as:
# lmeth1,lunmeth1,rmeth1,runmeth1,lmeth2,lunmeth2,rmeth2,runmeth2 etc
def covReplicateMethCounts(data,covs,compare):
    if len(compare)==2:
        covleft, covright = compare
        
        lreplicates = []
        rreplicates = []
    
        for cov, (meth,unmeth) in izip(covs,grouper(2,data)): # extract in pairs of meth/unmeth
            if cov == covleft:
                lreplicates.append((meth,unmeth))
            elif cov == covright:
                rreplicates.append((meth,unmeth))

        #return [(lmeth,lunmeth,rmeth,runmeth),(),(),...]    
        return [(a,b,c,d) for (a,b),(c,d) in izip(lreplicates,rreplicates)]
    else:
        # assume that each pair of samples (i.e. 4 columns) in the input is one replicate
        # in this case our input file is in the right order and we just need to split it up.
        return grouper(4,data)


# total G: at least some of the data is different
# 1 degree of freedom per replicate (i.e. 3 in total), add G's from each seperate test
def totalG(data,covs,compare):
    dof = 0
    totalg = 0
    for lmeth,lunmeth,rmeth,runmeth in covReplicateMethCounts(data,covs,compare): # extract in quads (2 pairs)
            totalg += gtest(lmeth,lunmeth,rmeth,runmeth)
            dof+=1 
    assert totalg >= 0
    return totalg, dof, GtoP(totalg,dof)

# pooled G: the data as a whole is different
# 1 degree of freedom, pool results from all replicates and do 1 2x2 table
def pooledG(data,covs,compare):    
    lmeth,lunmeth,rmeth,runmeth = pooledMethCounts(data,covs,compare)
    
    dof = 1
    
    pooledg = gtest(lmeth,lunmeth,rmeth,runmeth)
    
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
        opts, args = getopt.getopt(sys.argv[1:], "", ["combined=","cov=","compare=","window=","nogtest"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    windowSize = 2000
    infile = None
    covs = None
    compare = None
    doGtest = True
    
    for o, a in opts:
        if o=="--combined":
            infile = a
            print "combined", a
        elif o=="--cov":
            covs = a.split(",")
            print covs
        elif o=="--compare":
            compare = a.split(",")
            print compare
            assert len(compare)==2
        elif o=="--window":
            windowSize = int(a)
            print "window size:", a
        elif o=="--nogtest":
            doGtest = False
            print "Not doing G test"
    
    assert infile != None
    
    # we either have compare and covs or neither
    assert compare == None or len(compare)==2
    if compare == None: assert covs == None
    if compare != None and len(compare)==2: assert covs != None
    
    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    
    # ---    
    # windows
    
    if compare != None:
        leftcov,rightcov = compare
        if doGtest:
            windows = csv.writer(open(infile+"."+leftcov+"-vs-"+rightcov+".windows.bedGraph","w"),delimiter='\t')
        percent = csv.writer(open(infile+"."+leftcov+"-vs-"+rightcov+".percent.windows.bedGraph","w"),delimiter='\t')
    else:
        if doGtest:
            windows = csv.writer(open(infile+".windows.bedGraph","w"),delimiter='\t')
        percent = csv.writer(open(infile+".percent.windows.bedGraph","w"),delimiter='\t')
    

    
    currentWindowPos = None
    currentChr = None
    
    cs = 0
    dataCount = None
    
    def writeWindow(chrm,pos,epos):
        global dataCount
        global cs
        global windows
        global percent
        
        lmethCount,lunmethCount,rmethCount,runmethCount = pooledMethCounts(dataCount,covs,compare)
        
        if (lmethCount + lunmethCount) == 0 or (rmethCount + runmethCount) == 0:
            print "Unprintable window:",chrm, pos, epos, lmethCount, lunmethCount, rmethCount, runmethCount
            return
        
        lPercentage = (float(lmethCount)/float(lmethCount+lunmethCount)) * 100.0
        rPercentage = (float(rmethCount)/float(rmethCount+runmethCount)) * 100.0
        
        percent.writerow([chrm,pos,epos+1,lPercentage,rPercentage])
        
        if doGtest:
            total,totaldof,totalp = totalG(dataCount,covs,compare)
            pooled,pooleddof,pooledp = pooledG(dataCount,covs,compare)
            hetero,heterodof,heterop = heterogeneityG(total,totaldof,pooled,pooleddof)
            
            windows.writerow([chrm,pos,epos+1,'1' if rPercentage > lPercentage else '-1',totalp,pooledp,heterop,cs,lmethCount,lunmethCount,rmethCount,runmethCount])
        
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
                
            # write a window if the next cpg is outside the current window
            elif coord - currentWindowPos > windowSize:
                writeWindow(currentChr,currentWindowPos,lastCoordCounted)
                currentChr = chrm
                currentWindowPos = coord
                
            # otherwise dont write any window (the current cpg is in the window)
        
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