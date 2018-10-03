'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv
import math

def parseMethLine(line):
    (chr,coord, meth, unmeth) = line
    return line
    #(readid,meth,chr,coord,methcall) = line
    #return (True if meth == '+' else False, chr, coord)

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["bismark=",
                                                      "methBedGraph=",
                                                      "percentageBedGraph=",
                                                      "coverageBedGraph="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    percentageBedGraph = None
    coverageBedGraph = None 
    
    for o, a in opts:
        if o=="--bismark":
            infile = a
            print "Bismark", a
        elif o=="--methBedGraph":
            methBedGraph = a
            print "Meth Bed Graph", a
        elif o=="--percentageBedGraph":
            percentageBedGraph = a
            print "Percentage Bed Graph", a
        elif o=="--coverageBedGraph":
            coverageBedGraph = a
            print "Coverage Bed Graph", a
    
    assert infile != None
    assert methBedGraph != None
    assert percentageBedGraph != None
    assert coverageBedGraph != None
    
    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    
    methBedGraphFile = csv.writer(open(methBedGraph,"w"),delimiter='\t')
    percentageBedGraphFile = csv.writer(open(percentageBedGraph,"w"),delimiter='\t')
    coverageBedGraphFile = csv.writer(open(coverageBedGraph,"w"),delimiter='\t')
    
    minReads = 10
    
    def outputC(chr, coord, methylated, unmethylated):
        
        methScore = methylated - unmethylated
        
        # log the methScore (note you can't log a 0 or a negative number so do some fixing for that
        if methScore != 0:
            if methScore > 0:
                methScore = math.log(methScore)
            else:
                methScore = -1*math.log(abs(methScore))
        
        methBedGraphFile.writerow([chr,coord,coord+1,methScore])
        
        # only do percentage if at least minReads number of reads
        if methylated + unmethylated > minReads:        
            percentageBedGraphFile.writerow([chr, coord, coord+1, (float(methylated)/float(methylated+unmethylated)) * 100.0])
        
        coverageBedGraphFile.writerow([chr, coord, coord+1, methylated+unmethylated])
    
    currentchr = None
    currentpos = None
    
    methylated = 0
    unmethylated = 0
    
    for line in methfile:
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
            outputC(currentchr,currentpos,methylated,unmethylated)
            
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
    outputC(currentchr,currentpos,methylated,unmethylated)
        