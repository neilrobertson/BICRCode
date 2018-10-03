'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv
from scipy.stats import binom_test
from itertools import izip_longest
from multiprocessing import Pool

def binomRow(row):
    chrm, coord, data = row

    outputRow = [chrm, coord]
    for rowextension in map(binomMethCall,data):
        outputRow.extend(rowextension)
    return outputRow

def binomMethCall(group):
        meth,unmeth,errorrate = group
        return [meth,
                unmeth,
                binom_test(meth,meth+unmeth,errorrate)]

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def unpackLine(chrm,coord,*data):
    return chrm.strip(), int(coord), [int(d) for d in data]

def extractLine(line):
    chrm,coord,data = unpackLine(*line)
    if not chrm.startswith("chr"):
        chrm = "chr"+chrm
    return (chrm,coord,data)
    
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["combined=","errorRates=","cores="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    cores = 4
    
    for o, a in opts:
        if o=="--combined":
            infile = a
            print "combined", a
        elif o=="--errorRates":
            errorRates = [float(x) for x in a.split(",")]
            print "Error Rates", errorRates
        elif o=="--cores":
            cores = int(a)
    
    assert infile != None
    
    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    
    binomOutput = csv.writer(open(infile+".binom","w"),delimiter='\t')

    p = Pool(cores)
    
    rows = []    
    for count,line in enumerate(methfile):
        
        if count % 100000 == 0:
            print count
            for row in p.map(binomRow,rows):
                binomOutput.writerow(row)
            rows = []
        
        chrm,coord,data = extractLine(line)
        
        row = [chrm,coord]
        
        # extract in pairs
        samples = [(meth,unmeth,errorRates[sampleno]) for sampleno,(meth,unmeth) in enumerate(grouper(2,data))]        
        row.append(samples)
        rows.append(row)
    
    print count
    for row in p.map(binomRow,rows):
        binomOutput.writerow(row)
        