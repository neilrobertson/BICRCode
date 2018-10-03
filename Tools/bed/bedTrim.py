'''
Created on 16 Jun 2011

@author: mcbryan
'''

import sys
from bed.treatment import SimpleBed
from genemapping.chrmEnds import ChromosomeEnds
import csv
import getopt

if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "b:a:o:", [])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    bedfile = None
    outputfile = None
    for o, a in opts:
        if o == "-a":
            assembly = a
        if o == "-b":
            bedfile = SimpleBed(a)
        if o == "-o":
            outputfile = csv.writer(open(a, "w"), delimiter='\t')

    ends = ChromosomeEnds(assembly)

    assert bedfile != None
    assert outputfile != None
    
    
    
    for (chr, start, stop) in bedfile:
        if int(start) > ends[chr]:
            continue
        else:
            outputfile.writerow([chr,start,str(min(ends[chr],int(stop)))])