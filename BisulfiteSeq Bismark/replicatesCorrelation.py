'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv
import math
import collections
from itertools import izip_longest
from scipy.stats.stats import pearsonr
import gc
gc.disable()
    
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

def getReplicates(data,quads=True):
    replicates = []
    
    if quads:
        for lmeth,lunmeth,rmeth,runmeth in grouper(4,data): # extract in quads (2 pairs)
            replicates.append((lmeth,lunmeth,rmeth,runmeth))
    else:
        for meth,unmeth in grouper(2,data): # extract in pairs (1 sample)
            replicates.append((meth,unmeth))
    
    return replicates

    
if __name__ == '__main__':   
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["combined=","a=","b=","side="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    a = None
    b = None
    side = None
    
    for opt, arg in opts:
        if opt=="--combined":
            infile = arg
            print "combined", arg
        if opt=="--a":
            a = int(arg)-1
            print "A:", arg
        if opt=="--b":
            b = int(arg)-1
            print "B:", arg
        if opt=="--side":
            side = arg
            print "Side:", arg
            assert side in ["left","right","absolute"]
    
    assert infile != None

    
    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    
    # ---
    
    minReadsInAtLeastOneSampleForPercentage = 10
    minReadsInBothSamplesForPercentage = 8
    
    bins = collections.defaultdict(list)

    lines = 0
    
    x = list()
    y = list()

    for line in methfile:
        
        lines+=1
        
        if lines % 1000000 == 0:
            print str(lines / 1000000) + " million"
        
        chrm,coord,data = extractLine(line)
        
        replicates = getReplicates(data, quads = False if side == "absolute" else True)
        
        # this should work for both absolute and left/right positioning
        ameth, aunmeth = replicates[a][2:4] if side == "right" else replicates[a][0:2]
        bmeth, bunmeth = replicates[b][2:4] if side == "right" else replicates[b][0:2]
        
        # only write the percentages if at least minReads number of reads
        if (ameth + aunmeth >= minReadsInAtLeastOneSampleForPercentage or bmeth + bunmeth >= minReadsInAtLeastOneSampleForPercentage) and (ameth + aunmeth >= minReadsInBothSamplesForPercentage and bmeth + bunmeth >= minReadsInBothSamplesForPercentage):
            
            aPercentage = (float(ameth)/float(ameth+aunmeth)) * 100.0
            bPercentage = (float(bmeth)/float(bmeth+bunmeth)) * 100.0
            
            bins[int(aPercentage)].append(bPercentage)
            
            x.append(aPercentage)
            y.append(bPercentage)
    
    print "Pearson:"+ str(pearsonr(x,y))
    
    for i in bins:
        avg = math.fsum(bins[i])/len(bins[i])
        print i, avg, avg-float(i)