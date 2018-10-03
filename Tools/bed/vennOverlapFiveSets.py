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
        opts, args = getopt.getopt(sys.argv[1:], "a:b:c:d:e:n", ["cores="])
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
        elif opt == "-d":
            d = arg
        elif opt == "-e":
            e = arg    
        elif opt == "-n":
            numberTrials = int(arg)
        elif opt == "--cores":
            cores = int(arg)
    
    #Gets the overlap values
    aa = numbBP(SimpleBed(a))
    bb = numbBP(SimpleBed(b))
    cc = numbBP(SimpleBed(c))
    dd = numbBP(SimpleBed(d))
    ee = numbBP(SimpleBed(e))
    
    ab,abBed = compareTwoFiles(a,b)
    ac,acBed = compareTwoFiles(a,c)
    ad,adBed = compareTwoFiles(a,d)
    ae,aeBed = compareTwoFiles(a,e)
    bc,bcBed = compareTwoFiles(b,c)
    bd,bdBed = compareTwoFiles(b,d)
    be,beBed = compareTwoFiles(b,e)
    cd,cdBed = compareTwoFiles(c,d)
    ce,ceBed = compareTwoFiles(c,e)
    de,deBed = compareTwoFiles(d,e)
    
    abc,abcBed= compareListAndFile(abBed,c)
    abd,abdBed= compareListAndFile(abBed,d)
    abe,abeBed= compareListAndFile(abBed,e)
    acd,acdBed= compareListAndFile(acBed,d)
    ace,aceBed= compareListAndFile(acBed,e)
    ade,adeBed= compareListAndFile(adBed,e)
    bcd,bcdBed= compareListAndFile(bcBed,d)
    bce,bceBed= compareListAndFile(bcBed,e)
    bde,bdeBed= compareListAndFile(bdBed,e)
    cde,cdeBed= compareListAndFile(cdBed,e)
    
    abcd,abcdBed= compareListAndFile(abcBed,d)
    abce,abceBed= compareListAndFile(abcBed,e)
    abde,abdeBed= compareListAndFile(abdBed,e)
    acde,acdeBed= compareListAndFile(acdBed,e)
    bcde,bcdeBed= compareListAndFile(bcdBed,e)
    
    abcde,_= compareListAndFile(abcdBed,e)
    
    print "Raw:"
    print str(aa)
    print str(bb)
    print str(cc)
    print str(dd)
    print str(ee)
    print str(ab)
    print str(ac)
    print str(ad)
    print str(ae)
    print str(bc)
    print str(bd)
    print str(be)
    print str(cd)
    print str(ce)
    print str(de)
    print str(abc)
    print str(abd)
    print str(abe)
    print str(acd)
    print str(ace)
    print str(ade)
    print str(bcd)
    print str(bce)
    print str(bde)
    print str(cde)
    print str(abcd)
    print str(abce)
    print str(abde)
    print str(acde)
    print str(bcde)
    print str(abcde)
    print 
    
    bcde = bcde - abcde
    acde = acde - abcde
    abce = abce - abcde
    abcd = abcd - abcde
    abde = abde - abcde
    cde = cde - bcde - acde
    bde = bde - bcde - abde
    bce = bce - bcde - abce
    bcd = bcd - bcde - abcd
    ade = ade - abde - acde 
    ace = ace - acde - abce
    acd = acd - acde - abcd
    abe = abe - abce - abde
    abd = abd - abcd - abde
    abc = abc - abce - abcd
    ab = ab - abe - abd - abc
    ac = ac - ace - acd - abc
    ad = ad - ade - acd - abd
    ae = ae - ade - ace - abe
    bc = bc - bce - bcd - abc
    bd = bd - bcd - abd - bde
    be = be - bde - bce - abe
    cd = cd - cde - bcd - acd
    ce = ce - cde - bce - ace
    de = de - bde - ade - cde
    aa = aa - ab - ac - ad - ae
    bb = bb - ab - bc - bd - be
    cc = cc - ac - bc - cd - ce
    dd = dd - ad - bd - cd - de
    ee = ee - ae - be - ce - de
    
    print "Results:"
    print str(aa)
    print str(bb)
    print str(cc)
    print str(dd)
    print str(ee)
    print str(ab)
    print str(ac)
    print str(ad)
    print str(ae)
    print str(bc)
    print str(bd)
    print str(be)
    print str(cd)
    print str(ce)
    print str(de)
    print str(abc)
    print str(abd)
    print str(abe)
    print str(acd)
    print str(ace)
    print str(ade)
    print str(bcd)
    print str(bce)
    print str(bde)
    print str(cde)
    print str(abcd)
    print str(abce)
    print str(abde)
    print str(acde)
    print str(bcde)
    print str(abcde)
    
    print "Vennerable:"
    print str(aa) + "," + str(bb) + "," + str(ab) + "," + str(cc) + "," + str(ac) + "," + str(bc) + "," + str(abc)
    
def concatenateFiles (self, path,tempfiles):
     
    f = open(path + "temp.txt", "w")
    for tempfile in tempfiles:
        f.write(tempfile.read())
          
    return path + "temp.txt"
     
    