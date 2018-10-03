#!/usr/bin/env python

# Maps a list of genes to their location!!
import os
import sys
import csv
import getopt

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], " ", ["genes=","geneLocations=","locationIDField=","outputFile="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) 
        sys.exit(2)

    genesListFilename = None
    geneLocationsFilename = None
    idField = 1
    outputFilename = None
    
    delimiter = "\t"
    
    for o, a in opts:
        if o == "--genes":
            genesListFilename = a
        elif o == "--geneLocations":
            geneLocationsFilename = a
        elif o == "--locationIDField":
            idField = int(a)
        elif o == "--outputFile":
            outputFilename = a
            
    assert genesListFilename != None
    assert geneLocationsFilename != None
    assert outputFilename != None
    
    genesList = []
    with open(genesListFilename, "r") as genesListFile:
        for row in genesListFile:
            genesList.append(row.split(delimiter)[0])
    genesList = [x.strip().upper() for x in genesList]
            
    with open(geneLocationsFilename, "r") as geneLocationsFile:
        with open(outputFilename, "w") as outputFile:
            for row in geneLocationsFile:
                mappingID = row.split(delimiter)[idField -1].strip().upper()
                if mappingID in genesList:
                    outputFile.write(row) 
            