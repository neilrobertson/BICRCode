#!/usr/bin/env python

import os
import sys

# add executing directory as path
sys.path.append(sys.path[0])
sys.path.append(os.path.expanduser("~/mount/repository/shared/baselib"))

import csv
import getopt
import numpy
import math

from bed.treatment import Bed as BedTreatment
from bed.treatment import SimpleBed

# treatment behaviours

if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "o:", [])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    # add executing directory as part of path
    sys.path.append(sys.path[0])
    
    for o, a in opts:
        if o == "-o":
            outputFile = a
        else:
            assert False,  "Unhandled option"
    
    if len(args) <= 1:
        assert False, "Only one list"
  
    class Interval():
        def __init__(self, start, stop):
            self.start = start
            self.stop = stop
        
        def __repr__(self):
            return str(self.start) + "-" + str(self.stop)
  
    def mergeIntervals(intervals):
        for i in range(len(intervals)):
            for j in range(i, len(intervals)):
                if i == j:
                    continue
                else:
                    if (intervals[j].start < intervals[i].stop and intervals[i].start < intervals[j].stop):
                        intervals[i].start = min(intervals[i].start,  intervals[j].start)
                        intervals[i].stop = max(intervals[i].stop,  intervals[j].stop)
                        del intervals[j]
                        mergeIntervals(intervals)
                        return
  
    allintervals = {}
    
    for arg in args:
        treatment = SimpleBed(arg)
        for (chr, start, stop) in treatment:
            interval = Interval(start, stop)
            if chr not in allintervals:
                allintervals[chr] = []
            allintervals[chr].append(interval)
    
    for chr in allintervals:
        mergeIntervals(allintervals[chr])

    treatments = {}
    
    for arg in args:
        treatments[arg] = BedTreatment(arg)
    
    prettyheader = ""
    for arg in args:
        prettyheader = prettyheader +"\t"+arg
        
    print "seq \t "+prettyheader
    
    for chr in allintervals:
        for interval in allintervals[chr]:
            valid = True
            hits = []
            for arg in args:
                newhits = treatments[arg].getOverlappingIntervals(chr, interval.start, interval.stop)
                if len(newhits) == 0:
                    valid = False
                    break
                else:
                    hits.extend(newhits)
            if valid:
                prettyhits = ""
                for hit in hits:
                    prettyhits = prettyhits + "\t"+str(hit)
                print chr+":"+str(interval.start)+"-"+str(interval.stop) +"\t", prettyhits
