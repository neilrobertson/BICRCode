#!/usr/bin/env python

import gc
import sys
import getopt
from bed.treatment import Bed as BedIntervalTree
from genemapping import Ensembl
from bed.treatment import SimpleBed

gc.disable()

try:
    opts, args = getopt.getopt(sys.argv[1:], "", [])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    sys.exit(2)


UPSTREAM_PROMOTOR_DIST = 5000
DOWNSTREAM_PROMOTOR_DIST = 1000


# probably want to change this to be exons rather than genes


###

# load data

genedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")

bp = 0
for gene in genedata:
    bp += genedata[gene].end -genedata[gene].start 

print bp

#genedata = Ensembl.GenesMapping(os.path.expanduser("~/mount/publicdata/hg18/ncbi36.1/genes-and-exons-human-NCBI36.1.csv"))

genes = Ensembl.ReverseGeneMapping(genedata)

genespluspromotor = Ensembl.ReverseGeneMapping(genedata, tssPadding = UPSTREAM_PROMOTOR_DIST)

genepromotors = Ensembl.ReversePromotorMapping(genedata, upstreamPadding = UPSTREAM_PROMOTOR_DIST, downstreamPadding = DOWNSTREAM_PROMOTOR_DIST)

exons = Ensembl.ReverseExonMapping(genedata)

class Interval():
    def __init__(self,start,end):
        self.start = start
        self.end = end

def mergeIntervals(intervals):
    sorted(intervals, key=lambda element: (element.start, element.end))
    
    mergedIntervals = []
    
    lastInterval = None
    for interval in intervals:
        if lastInterval != None:
            if (lastInterval.start < interval.end and interval.start < lastInterval.end):
                # overlapping, extend current interval
                interval.start = min(lastInterval.start,  interval.start)
                interval.end = max(lastInterval.end,  interval.end)
            else:
                # not overlapping, add last Interval to mergedIntervals
                mergedIntervals.append(lastInterval)
        # last interval becomes current interval
        # (either extended from the overlap, or new from the list when not overlapping)
        lastInterval = interval
    
    if lastInterval != None:
        mergedIntervals.append(lastInterval)
    
    return mergedIntervals
                
def intersectingIntervals(aIntervals,bIntervals):
    intersections = []
    for ai in mergeIntervals(aIntervals):
        for bi in mergeIntervals(bIntervals):
            if (ai.start < bi.end and bi.start < ai.end):
                intersections.append(Interval(max(ai.start,bi.start), min(ai.end,bi.end)))
    return intersections


def overlappingBP(start,stop,intervals):
    overlapBP = 0
    for interval in intervals:
        # how much of interval is overlapped with start, stop
        overlapBP += max(min(interval.end,stop) - max(interval.start,start),0)
    return overlapBP

for infile in args:

    print infile

    intervals = SimpleBed(infile)
    
    print "Starting..."
    
    for row in intervals:
        
        chr, start, stop = row
        
        # 0 length regions
        if start == stop:
            continue
        
        size = stop - start
        
        # genes
        ingenes = mergeIntervals(genes.getIntervalsInRange(chr, start, stop))
        geneoverlap = overlappingBP(start,stop,ingenes)
        
        # genes with promotors
        ingenespluspromotors = mergeIntervals(genespluspromotor.getIntervalsInRange(chr, start, stop))
        genepluspromotoroverlap = overlappingBP(start,stop,ingenespluspromotors)
        
        # promotors
        inpromotors = mergeIntervals(genepromotors.getIntervalsInRange(chr, start, stop))
        promotoroverlap = overlappingBP(start,stop,inpromotors)
        
        inexons = mergeIntervals(exons.getIntervalsInRange(chr, start, stop))
        exonoverlap = overlappingBP(start,stop,inexons)
        
        intronoverlap = geneoverlap - exonoverlap
        
        # both gene and promotor
        bothIntersection = intersectingIntervals(ingenes,inpromotors)
        intersectionOverlap = overlappingBP(start,stop,bothIntersection)
        
        intergenic = size - geneoverlap - promotoroverlap + intersectionOverlap
        
        print chr, start, stop, size, geneoverlap, promotoroverlap, exonoverlap, intergenic, intronoverlap
    
exit(1)
