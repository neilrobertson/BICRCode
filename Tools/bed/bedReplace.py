#!/usr/bin/env python

import os
import sys
import csv
import getopt
import numpy
import math

from bed.treatment import Bed as BedTreatment

# assume vstep
debug = False

# treatment behaviours

if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "da:b:o:", [])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    # add executing directory as part of path
    sys.path.append(sys.path[0])
    
    for o, a in opts:
        if (o=="-d"):    
            debug = True

    for o, a in opts:
        if (o=="-d"):
            pass # we dealt with this already
        elif o == "-b":
            treatmentFileName = a
            if debug:           
                replaceTreatment = BedTreatment(treatmentFileName, 50000)
            else:
                replaceTreatment = BedTreatment(treatmentFileName)
        elif o == "-a":
            baseTreatmentFileName = a
        elif o == "-o":
            outputFile = a
        else:
            assert False,  "Unhandled option"
  
    baseTreatment = csv.reader(open(baseTreatmentFileName, "r"), delimiter='\t')
    outputFile = csv.writer(open(outputFile, "w"), delimiter='\t')
    
    for row in baseTreatment:
        if len(row)==0 or row[0].startswith("#"):
            continue
        chr = row[0]
        start = int(row[1])
        stop = int(row[2])
        intervals = replaceTreatment.getOverlappingIntervals(chr, start, stop)
        
        for interval in intervals:
            outputFile.writerow((chr, interval.start, interval.end))
