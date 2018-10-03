#!/usr/bin/env python

import os
import sys
import getopt
import csv
import math

def getChromEnds(filename):
    chromends = {}
    reader = csv.reader(open(filename), delimiter="\t")
    for row in reader:
        chromname = row[0]
        if not chromname.startswith("chr"):
            chromname = "chr" + chromname
        chromends[chromname] = int(row[1])
    return chromends

if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:i:o:e:", [])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    vstepWidth = 200 # default
    
    for o, a in opts:
        # need a vstep width
        if o=="-s":
            vstepWidth = int(a)
        elif o == "-i":
            inputFilename = a
        elif o == "-o":
            outputFilename = a
        elif o == "-e":
            chrmEnds = getChromEnds(a)
    
    # first work out the average
    inputFile = csv.reader(open(inputFilename), delimiter="\t")
    
    total = 0
    seenValues = 0
    seenChrms = set()
    
    for row in inputFile:
        (chr,  pos, value) = row
        total += int(value)
        seenValues += 1
        seenChrms.add(chr)
    
    # work out how many values we should have for a full vstep step file (only consider chrms we found something on)
    bases = 0
    for chrm in seenChrms:
        bases+= chrmEnds[chrm]
    expectedValues = bases/vstepWidth if bases%vstepWidth == 0 else (bases/vstepWidth) +1  # account for rounding
    
    # and we have an average - yay
    mean = float(total) / float(expectedValues)
    
    print "Mean:"+ str(mean)
    
    inputFile = csv.reader(open(inputFilename), delimiter="\t")
    
    squareddifftotal = 0
    seenChrms = set()
    
    for row in inputFile:
        (chr,  pos, value) = row
        
        difference = int(value) - mean

        squareddiff = difference * difference
        
        squareddifftotal += squareddiff

        seenChrms.add(chr)
    
    # now adjust the square diff total for any zeros we didnt see
    numberMissingValues = expectedValues - seenValues
    diffMissing = 0.0 - mean
    squareddiffMissing = diffMissing * diffMissing
    
    squareddifftotal = squareddifftotal + (squareddiffMissing * numberMissingValues)
    
    # now we have a stddev yay
    stddev = math.sqrt(float(squareddifftotal) / float(expectedValues))
    
    print "Std dev:" + str(stddev)
    
    # which means we can finally work out a zscore
    
    inputFile = csv.reader(open(inputFilename), delimiter="\t")
    outputFile = csv.writer(open(outputFilename, "w"), delimiter="\t")
    
    for row in inputFile:
        (chr,  pos, value) = row
        
        value = int(value)
        
        zscore = (float(value) - mean) / stddev
        
        outputFile.writerow([chr, pos, zscore])
