'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv
from scipy.stats import binom_test

def parseMethLine(line):
    (chrm,coord,lmeth,lunmeth) = line
    chrm = chrm.strip()
    if not chrm.startswith("chr"):
        chrm = "chr"+chrm
    # coord should be in base 0 already, just need it as an int
    coord = int(coord)
    meth=int(lmeth)
    unmeth=int(lunmeth)
    return (chrm,coord,meth,unmeth)

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["file=","errorRate="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    leftsuffix = None
    rightsuffix = None 
    
    for o, a in opts:
        if o=="--file":
            infile = a
            print "combined", a
        elif o=="--errorRate":
            errorRate = float(a)
            print "Error Rate", a

    assert infile != None
    
    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    binomOutput = csv.writer(open(infile+".binom","w"),delimiter='\t')
    
    for line in methfile:
        (chrm,coord,meth,unmeth) = parseMethLine(line)
        binomOutput.writerow([chrm,coord,
                              meth,unmeth,
                              binom_test(meth,meth+unmeth,errorRate)])