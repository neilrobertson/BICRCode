'''
Created on 18 Mar 2011

@author: mcbryan
'''

import sys
import getopt
import csv
import collections
from csvfile.indexedcsv import IndexedCSV,ColumnIndex
from genemapping import Ensembl

if __name__ == '__main__':
    
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", [
                                                      # command args go here
                                                      "gene-expression-difference=",
                                                      "read-group-tracking=",
                                                      "output="
                                                      ])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
    for o, a in opts:
        if o=="--gene-expression-difference":
            differences = IndexedCSV(a, key="test_id")
        elif o=="--read-group-tracking":
            readgrouptracking = IndexedCSV(a, useRowNumAsKey=True)
            xlocindex = ColumnIndex(readgrouptracking,"tracking_id")
        elif o=="--output":
            outputFile = a
        else:
            print "Unknown parameter: "+o+" "+a
            sys.exit(2)
    
    # get a suitable header
    replicatenames = []
    for readtrackingrow in xlocindex["XLOC_000001"]:
        condition = readgrouptracking[readtrackingrow]["condition"]
        replicate = readgrouptracking[readtrackingrow]["replicate"]
        name = condition+"."+replicate
        replicatenames.append(name)
    print replicatenames
    
    header = differences.keys[:]
    header.extend(replicatenames)
    print header
    
    output = csv.writer(open(outputFile,"w"),delimiter="\t")
    output.writerow(header)
    
    for xloc in differences:
        trackings = {}
        for readtrackingrow in xlocindex[xloc]:
            condition = readgrouptracking[readtrackingrow]["condition"]
            replicate = readgrouptracking[readtrackingrow]["replicate"]
            name = condition+"."+replicate
            fpkm = readgrouptracking[readtrackingrow]["FPKM"]
            trackings[name] = fpkm
        
        row = differences.rows[xloc]
        row.extend([trackings[name] for name in replicatenames])
        output.writerow(row)
