#!/usr/bin/env python

import sys
import getopt
import csv
from csvfile.indexedcsv import IndexedCSV

# Joins two CSV files based on a column
# Only joins one line per item
# Only joins keys which are in both files

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "1:2:o:", ["key1=","key2="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    location1 = None
    key1 = 0
    
    location2 = None
    key2 = 0
    
    outlocation = None
    
    for o, a in opts:
        if o=="-1":
            location1 = a
        elif o=="-2":
            location2 = a 
        elif o=="-o":
            outlocation = a
        elif o=="--key1":
            key1 = int(a)
        elif o=="--key2":
            key2 = int(a)
    
    one = IndexedCSV(location1,keyPos=key1)
    
    two = IndexedCSV(location2,keyPos=key2)
    
    with open(outlocation,"w") as outfile:
        csvout = csv.writer(outfile,delimiter="\t")
        
        header = list(one.keys)
        header.extend(two.keys)
        header.insert(0,"Key")
        
        csvout.writerow(header)        

        for item in one:
            
            row = [item]
            for key in one.keys:
                row.append(one[item][key])
            for key in two.keys:
                row.append(two[item][key])
            csvout.writerow(row)

