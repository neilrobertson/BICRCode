'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv
import math
import collections

def reverseOrder(line):
    (chrm,coord,lmeth,lunmeth,rmeth,runmeth) = line
    return (chrm,coord,rmeth,runmeth,lmeth,lunmeth)

def parseMethLine(line,reverse):
    if reverse:
        line = reverseOrder(line)
    
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

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["combined=","reverseOrder"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None

    reverse = False
    
    for o, a in opts:
        if o=="--combined":
            infile = a
            print "combined", a
        elif o=="--reverseOrder":
            reverse = True
    
    assert infile != None

    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    
    # ---
    
    minReadsInAtLeastOneSampleForPercentage = 10
    minReadsInBothSamplesForPercentage = 3
    
    bins = collections.defaultdict(list)

    lines = 0

    for line in methfile:
        
        (chrm,coord,lmeth,lunmeth,rmeth,runmeth) = parseMethLine(line,reverse)
        
        lines+=1
        
        if lines % 1000000 == 0:
            print lines / 1000000
        
        # only write the percentages if at least minReads number of reads
        if (lmeth + lunmeth >= minReadsInAtLeastOneSampleForPercentage or rmeth + runmeth >= minReadsInAtLeastOneSampleForPercentage) and (lmeth + lunmeth >= minReadsInBothSamplesForPercentage and rmeth + runmeth >= minReadsInBothSamplesForPercentage):
            
            lPercentage = (float(lmeth)/float(lmeth+lunmeth)) * 100.0
            rPercentage = (float(rmeth)/float(rmeth+runmeth)) * 100.0
            
            bins[int(lPercentage)].append(rPercentage)
    
    for i in bins:
        avg = math.fsum(bins[i])/len(bins[i])
        print i, avg, avg-float(i)