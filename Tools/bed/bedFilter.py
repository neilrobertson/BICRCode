#!/usr/bin/env python

import os
import sys

# add executing directory as path
sys.path.append(sys.path[0])
sys.path.append(os.path.expanduser("~/mount/repository/shared/baselib"))

import csv
import getopt

from bed.treatment import Bed as BedTreatment

# assume vstep
debug = False

# treatment behaviours

def usage():
    print "usage: bedFilter.py [SWITCHES]"
    print "all switches are mandatory except -p"
    print "	-a <filename.bed> file to filter"
    print "	-b <filename.bed> file to filter by (subtract)"
    print "	-p <padvalue> (int: default 0) extend the size of the entries in the subtract file in both directions"
    print "	-o <filename> output filename"

if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:b:o:p:", [])
    except getopt.GetoptError, err:
        # print help information and exit:
#        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    # add executing directory as part of path
    sys.path.append(sys.path[0])
    
    treatmentFileName = None
    baseTreatmentFileName = None
    outputFile = None
    pad = 0

    for o, a in opts:
        if o == "-b":
            treatmentFileName = a
            subtractTreatment = BedTreatment(treatmentFileName)
        elif o == "-a":
            baseTreatmentFileName = a
        elif o == "-o":
            outputFile = a
        elif o == "-p":
            pad = int(a)
        else:
            assert False,  "Unhandled option"
  
    if treatmentFileName is None or baseTreatmentFileName is None or outputFile is None:
        usage()
        sys.exit(2)

    baseTreatment = csv.reader(open(baseTreatmentFileName, "r"), delimiter='\t')
    outputFile = csv.writer(open(outputFile, "w"), delimiter='\t')
    
    for row in baseTreatment:
        if len(row)==0 or row[0].startswith("#"):
            continue
        chr = row[0]
        start = int(row[1]) - pad
        stop = int(row[2]) + pad
        intervals = subtractTreatment.getOverlappingIntervals(chr, start, stop)
        
        # just drop intervals that have no overlaps entirely for now
        if len(intervals)==0:
            continue
        
        outputFile.writerow(row)
