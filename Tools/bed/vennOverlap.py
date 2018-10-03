#!/usr/bin/env python

import sys
import getopt
import gc

gc.disable()

from bed.treatment import SimpleBed
from bed.treatment import Bed as BedIntervalTree

from datastructures.genomeintervaltree import GenomeIntervalTree


def determineOverlap(regions,treatment):
   
    overlappingIntervals = []    
    overlappingBP = 0
    
    for row in regions:
        (chr, start, stop) = row
        intervals = treatment.getIntervalsInRange(chr, start, stop)
            
        for interval in intervals:
            # how much of interval is overlapped with chr, start, stop
            overlappingBP += min(interval.end,stop) - max(interval.start,start)
            
            overlappingIntervals.append((chr,max(interval.start,start),min(interval.end,stop)))
            
    return overlappingBP, overlappingIntervals

def numbBP(regions):
    bp = 0
    for chr,start,stop in regions:
        bp += stop-start
    return bp


def compareTwoFiles(a,b):
    
    aBed = SimpleBed(a)
    bTree = BedIntervalTree(b)

    return determineOverlap(aBed,bTree)
    
def compareListAndFile(ab,c):
    cTree = BedIntervalTree(c)
    
    return determineOverlap(ab,cTree)
    
# treatment behaviours

if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:b:c:n", ["cores="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    # add executing directory as part of path
    sys.path.append(sys.path[0])

    numberTrials = 1000
    cores = 1

    for opt, arg in opts:
        if opt == "-a":
            a = arg
        elif opt == "-b":
            b = arg
        elif opt == "-c":
            c = arg    
        elif opt == "-n":
            numberTrials = int(arg)
        elif opt == "--cores":
            cores = int(arg)
    
    #Gets the overlap values
    aa = numbBP(SimpleBed(a))
    
    bb = numbBP(SimpleBed(b))
    
    cc = numbBP(SimpleBed(c))
    
    ab,abBed = compareTwoFiles(a,b)
    
    ac,_ = compareTwoFiles(a,c)
    
    bc,_ = compareTwoFiles(b,c)
    
    abc,_= compareListAndFile(abBed,c)
    
    print "Raw:"
    print "a\t" + str(aa) 
    print "b\t" + str(bb) 
    print "c\t" + str(cc) 
    print "ab\t" + str(ab) 
    print "ac\t" + str(ac) 
    print "bc\t" + str(bc) 
    print "abc\t" + str(abc) + "\n"
    
    ab = ab - abc
    ac = ac - abc
    bc = bc - abc
    aa = aa - ab - ac - abc
    bb = bb - ab - bc - abc
    cc = cc - ac - bc - abc
    
    print "Results:"
    print "a\t" + str(aa)
    print "b\t" + str(bb)
    print "ab\t" + str(ab)
    print "c\t" + str(cc)
    print "ac\t" + str(ac)
    print "bc\t" + str(bc)
    print "abc\t" + str(abc)
    
    print "Vennerable:"
    print str(aa) + "," + str(bb) + "," + str(ab) + "," + str(cc) + "," + str(ac) + "," + str(bc) + "," + str(abc)
    
    
    
    