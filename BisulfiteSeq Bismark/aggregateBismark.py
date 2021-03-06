'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv

def parseMethLine(line):
    (readid,meth,chr,coord,methcall) = line
    return (True if meth == '+' else False, chr, coord)

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   "",
                                   ["bismark=", "aggregate="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    aggregatefile= None
    
    for o, a in opts:
        if o=="--bismark":
            infile = a
            print "Bismark", a
        elif o=="--aggregate":
            aggregatefile = a
            print "Aggregate", a
    
    assert infile != None
    assert aggregatefile != None
    
    bismark = csv.reader(open(infile, "r"), delimiter="\t")
    
    aggregate = csv.writer(open(aggregatefile,"w"),delimiter='\t')      
    
    currentchr = None
    currentpos = None
    
    methylated = 0
    unmethylated = 0
    
    for line in bismark:
        (methcall, chr, coord) = parseMethLine(line)
        
        if not chr.startswith("chr"):
            chr = "chr"+chr
        
        coord = int(coord) -1 # meth calls in 1 base, ucsc / bedgraph in 0 base
        
        # init
        if currentpos == None or chr == None:
            currentchr = chr
            currentpos = coord
        elif currentpos != coord or currentchr != chr:
            # new one
            # output our current status
            aggregate.writerow([currentchr,currentpos,methylated,unmethylated])
            
            # reset to fresh
            currentchr = chr
            currentpos = coord
            methylated = 0
            unmethylated = 0
        
        # work at current location
        if methcall:
            methylated += 1
        else:
            unmethylated += 1
    
    # output last line
    aggregate.writerow([currentchr,currentpos,methylated,unmethylated])
    