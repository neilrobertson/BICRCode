'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv
from sequence.genome import Genome

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["bedgraph="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    aggregatefile= None
    
    for o, a in opts:
        if o=="--bedgraph":
            infile = a
            print "Bedgraph", a
    
    assert infile != None
    
    build="hg18"
    
    genome = Genome(genomeBuild = build)

    def isForwardStrand(line):
        (chr,start,stop,value) = line
        base = genome.getSequence(chr,int(start),int(stop)).upper()
        assert base in ["C","G"]
        return True if base == "C" else False

    def isReverseStrand(line):
        return not forwardStrand(line)
    
    bedgraph = csv.reader(open(infile, "r"), delimiter="\t") 
    forwardStrand = csv.writer(open(infile+"-forward.bedGraph","w"),delimiter='\t')
    reverseStrand = csv.writer(open(infile+"-reverse.bedGraph","w"),delimiter='\t')
    
    for line in bedgraph:        
        if isForwardStrand(line):
            # forward strand here
            forwardStrand.writerow(line)
        else:
            # reverse strand here
            reverseStrand.writerow(line)