#!/usr/bin/env python

"""
The plan here was to use the hypergeometric distribution to measure significance in the overlap.
However, that seems to be a bust as it's highly sensitive to the number of samples.

In this case the number of samples is huge (size of genome) so virtually everything is significant even if
only a fc of 1.000001 etc

If we return to this consider the methods in : http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3475982/#R44
"""

import sys
import getopt

from bed.treatment import SimpleBed
from bed.treatment import Bed as BedIntervalTree

from genemapping.gaps import ChromosomeGaps
from genemapping.chrmEnds import ChromosomeEnds
from genemapping.chrList import ChrList

from datastructures.genomeintervaltree import GenomeIntervalTree

import scipy.stats

def buildIntervalTree(self):
    intervalTree = GenomeIntervalTree()
    for chrm,start,stop in self:
        intervalTree.insertInterval(chrm, start, stop, None)
    return intervalTree

def determineOverlap(a,b):
    bTree = buildIntervalTree(b)

    totalOverlappingBP = 0
    
    for chrm, start, stop in a:
        intervals = bTree.getIntervalsInRange(chrm, start, stop)
            
        for interval in intervals:
            # how much of interval is overlapped with chr, start, stop
            totalOverlappingBP += min(interval.end,stop) - max(interval.start,start)
        
    return totalOverlappingBP

def numbBP(regions):
    bp = 0
    for chrm,start,stop in regions:
        bp += stop-start
    return bp

"""
   Calculates the pvalue using the hypergeometric distribution
   c: bp in common
   t: total bp
   a: bp in A
   b: bp in B
"""
def pValue(c,a,b,t):

    print c,a,b,t
    t = t - a - b + c
    a = a - c
    b = b - c


    if sum([c,t,a,b]) >= 100000:
        import rpy2.rpy_classic as rpy
        import rpy2.robjects as robjects
        
        m = robjects.r("matrix(c("+str(c)+","+str(a)+","+str(b)+","+str(t)+"),byrow=TRUE,nrow=2)")
        rpy.set_default_mode(rpy.BASIC_CONVERSION)

        val = rpy.r.chisq_test(m)
        return val
    else:
        return scipy.stats.fisher_exact([[c,a],[b,t]])

def compareTwo(a,b,effectiveGenomeSize,validChrs):
    
    print "Compare 2"
    aBed = SimpleBed(a,validChrs=validChrs)
    bBed = SimpleBed(b,validChrs=validChrs)

    print "Determine overlap"
    overlappingbp = determineOverlap(aBed,bBed)
    bpA = numbBP(aBed)
    bpB = numbBP(bBed)


    print "A percent overlapped:"+str(100.0*float(overlappingbp)/float(bpA)) + ", Expected:"+ str(100.0*float(bpB)/float(effectiveGenomeSize))
    print "B percent overlapped:"+str(100.0*float(overlappingbp)/float(bpB)) + ", Expected:"+ str(100.0*float(bpA)/float(effectiveGenomeSize))
    
    print "Pvalue" + pValue(overlappingbp,bpA,bpB,effectiveGenomeSize)
    

if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:b:d", ["assembly="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    # add executing directory as part of path
    sys.path.append(sys.path[0])

    direction = True

    assembly = "hg18"

    for opt, arg in opts:
        if opt == "-a":
            a = arg
        elif opt == "-b":
            b = arg
        elif opt == "--assembly":
            assembly = arg
    
    gaps = ChromosomeGaps(assembly)
    ends = ChromosomeEnds(assembly)
    
    effectiveGenomeSize = 0
    for chr in ChrList(assembly):
        effectiveGenomeSize += ends[chr] - gaps.chrmgaps[chr]
    
    print "---"
    print a.split("/")[-1] + " which have " + b.split("/")[-1]
    print "---"
    
    compareTwo(a,b,effectiveGenomeSize,ChrList(assembly))
