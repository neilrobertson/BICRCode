'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv
from itertools import izip_longest

class CommentException(Exception):
    pass

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def unpackLine(chrm,coord,*data):
    return chrm.strip(), int(coord), [int(d) for d in data]

def extractLine(line):
    if line[0].startswith("#"):
        raise CommentException()
    chrm,coord,data = unpackLine(*line)

    if not chrm.startswith("chr"):
        chrm = "chr"+chrm
    return (chrm,coord,data)
    
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["combined=","suffixes="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    suffixes = None
    rightsuffix = None
    
    for o, a in opts:
        if o=="--combined":
            infile = a
            print "combined", a
        elif o=="--suffixes":
            suffixes = a.split(",")
            print "Suffix", suffixes

    
    assert infile != None
    assert suffixes != None
    
    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    
    files = {}
    
    for suffix in suffixes:
        percentage = csv.writer(open(infile+"-"+suffix+".percent.bedGraph","w"),delimiter='\t')
        coverage = csv.writer(open(infile+"-"+suffix+".coverage.bedGraph","w"),delimiter='\t')
        files[suffix]=(percentage,coverage)
    
    # ---
    # for each
    # coverage
    # percentage
    
    minReadsForPercentage = 5
    minSamplesWithCoverage = 1

    for line in methfile:
        try:
            (chrm,coord,data) = extractLine(line)
        except CommentException:
            continue
        
        percentages = []
        coverages = []
        
        for meth,unmeth in grouper(2,data):
            coverage = meth + unmeth
            if coverage == 0:
                percentage = 0
            else:           
                percentage = (float(meth)/float(meth+unmeth)) * 100.0
            percentages.append(percentage)
            coverages.append(coverage)
    
        satisfycoveragecriteria = 0
        for coverage in coverages:
            if coverage >= minReadsForPercentage:
                satisfycoveragecriteria += 1
        if satisfycoveragecriteria < minSamplesWithCoverage:
            continue

        for index,(percentage,coverage) in enumerate(izip_longest(percentages,coverages)):
            
            percentfile,coveragefile = files[suffixes[index]]
            percentfile.writerow([chrm,coord,coord+1,percentage])
            coveragefile.writerow([chrm,coord,coord+1,coverage])
            