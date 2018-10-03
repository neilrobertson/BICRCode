#!/usr/bin/env python

# Grep output just for BP overlaps in both directions:
# cat tmp | grep -E "^(BP|Overlap BP)" | sed "s/.*(//" | sed 's/)//' | paste -d"," - - - -

import sys
import getopt
import gc

gc.disable()

from bed.treatment import SimpleBed
from bed.treatment import Bed as BedIntervalTree

from datastructures.genomeintervaltree import GenomeIntervalTree

from genemapping.gaps import ChromosomeGaps
from genemapping.chrmEnds import ChromosomeEnds
from genemapping.chrList import ChrList

from multiprocessing import Pool

import random
import collections

def randomlyAllocateRegions(regions):
    
    treatment = GenomeIntervalTree()

    for row in regions:
        
        (chr, start, stop) = row
        length = stop - start
        
        # randomly select chromosome
        #chr = randomChromosome()
        
        # find a random start/stop position that doesnt overlap with a gap
        suitable = False
        while not suitable:
            start = random.random() * ends[chr]
            stop = start + length
            
            if len(gaps.getIntervalsInRange(chr,start,stop)) == 0:
                suitable = True
        
        # insert into the the treatment
        treatment.insertInterval(chr, start, stop, None)
        
    return treatment

def determineOverlap(regions,treatment):
    overlapping = 0
    nonoverlapping = 0
    totalOverlappingBP = 0
    percentOverlap = 0
    
    overlappingLengths = []
    nonOverlappingLengths= []
    
    for row in regions:
        (chr, start, stop) = row
        intervals = treatment.getIntervalsInRange(chr, start, stop)
        
        overlappingBP = 0
        
        if len(intervals)>0:
            overlapping += 1
            overlappingLengths.append(stop-start)
        else:
            nonoverlapping += 1
            nonOverlappingLengths.append(stop-start)
            
        for interval in intervals:
            # how much of interval is overlapped with chr, start, stop
            overlappingBP += min(interval.end,stop) - max(interval.start,start)
        
        totalOverlappingBP += overlappingBP
        
        if overlappingBP / float(stop-start) >= 0.9:
            percentOverlap += 1
            
    return overlapping,totalOverlappingBP,percentOverlap,0 if overlapping == 0 else sum(overlappingLengths)/float(overlapping),0 if nonoverlapping == 0 else sum(nonOverlappingLengths)/float(nonoverlapping)

def numbBP(regions):
    bp = 0
    for chr,start,stop in regions:
        bp += stop-start
    return bp

def doOverlap(p):
    aBed,bBed = p
    randomRegions = randomlyAllocateRegions(bBed)
    return determineOverlap(aBed,randomRegions)

def pValue(actual,randomTrials):
    gt = 0
    lt = 0
    for trial in randomTrials:
        if trial >= actual:
            gt += 1
        elif trial <= actual:
            lt += 1
    return gt/float(len(randomTrials)),lt/float(len(randomTrials))


def compareTwo(a,b,numbRandomRegions,cores, validChrs):
    
    randomOverlap = []
    randomBPOverlap = []
    randomPercentOverlap = []
    
    aBed = SimpleBed(a,validChrs=validChrs)
    bBed = SimpleBed(b,validChrs=validChrs)
    
    params = (aBed,bBed)
    pool = Pool(processes = cores)
    result = pool.map(doOverlap, [params]*numbRandomRegions)
    
    for r,rbp,rpercent,_,_ in result:
        randomOverlap.append(r)
        randomBPOverlap.append(rbp)
        randomPercentOverlap.append(rpercent)
        

    bTree = BedIntervalTree(b)

    overlapping, overlappingbp, overlappingPercent,avgOverlappingLength,avgNonOverlappingLength = determineOverlap(aBed,bTree)
    
    #########
    
    random = sum(randomOverlap)/float(len(randomOverlap))
    
    print "Average Overlapping Length:"+str(avgOverlappingLength)
    print "Average NonOverlapping Length:"+str(avgNonOverlappingLength)
    
    print
    
    print "Overlap 1bp: " + str(overlapping) + " / " + str(len(aBed)) + " ("+str(100.0*overlapping/float(len(aBed)))+"%)"
    print "Overlap 1bp Random: " + str(random) + " / " + str(len(aBed)) + " ("+str(100.0*random/float(len(aBed)))+"%)"
    print "FC 1bp: "+str(overlapping) + " / " + str(random) + " = " + str(overlapping/float(random))
    
    #print "Max FC 1bp:"+ str(len(aBed)) + " / " + str(random) + " = " + str(len(aBed)/float(random))
    
    print "Pvalues (corr/anticorr):" + str(pValue(overlapping,randomOverlap))
    
    ##########
    
    print
    
    ##########
    
    randombp = sum(randomBPOverlap)/float(len(randomBPOverlap))
    
    print "BP: "+ str(overlappingbp) + " / " + str(numbBP(aBed)) + " ("+str(100.0*overlappingbp/float(numbBP(aBed)))+"%)"
    print "Overlap BP Random: " + str(randombp) + " / " + str(numbBP(aBed)) + " ("+str(100.0*randombp/float(numbBP(aBed)))+"%)"
    print "FC BP: "+str(overlappingbp) + " / " + str(randombp) + " = " + str(overlappingbp/float(randombp))

    #print "Max FC BP:"+ str(numbBP(aBed)) + " / " + str(randombp) + " = " + str(numbBP(aBed)/float(randombp))

    print "Pvalues (corr/anticorr):" + str(pValue(overlappingbp,randomBPOverlap))
    
    

# treatment behaviours

if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:b:n:d", ["cores=","assembly="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    # add executing directory as part of path
    sys.path.append(sys.path[0])

    numberTrials = 100
    cores = 16
    direction = True

    assembly = "hg18"

    for opt, arg in opts:
        if opt == "-a":
            a = arg
        elif opt == "-b":
            b = arg
        elif opt == "-n":
            numberTrials = int(arg)
        elif opt == "--cores":
            cores = int(arg)
        elif opt == "-d":
            direction = False
        elif opt == "--assembly":
            assembly = arg
    
    gaps = ChromosomeGaps(assembly)
    ends = ChromosomeEnds(assembly)
    
    from datastructures.wrg import WeightedRandomGenerator
    
    #weights = collections.defaultdict(int)
    
    effectiveGenomeSize = 0
    for chr in ChrList(assembly):
        effectiveGenomeSize += ends[chr] - gaps.chrmgaps[chr]
        # below is needed for weighted chromosome selection
        #weights[chr] = ends[chr] - gaps.chrmgaps[chr]
        #assert weights[chr] > 0
    print "Effective Genome Size:" + str(effectiveGenomeSize)
    
    # optional - removed for now.  allows regions to move between chromosomes
    #randomChromosome = WeightedRandomGenerator(weights)
    

    
    print "---"
    print a.split("/")[-1] + " which have " + b.split("/")[-1]
    print "---"
    
    compareTwo(a,b,numberTrials,cores,ChrList(assembly))
    
    print
    
    if direction:
        print "---"
        print b.split("/")[-1] + " which have " + a.split("/")[-1]
        print "---"
        
        compareTwo(b,a,numberTrials,cores,ChrList(assembly))