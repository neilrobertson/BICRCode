'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv
import math
# https://github.com/brentp/fishers_exact_test 
from fisher import cfisher as fisher_exact
from scipy.stats import binom_test

def parseMethLine(line):
    (chrm,coord,lmeth,lunmeth,rmeth,runmeth) = line
    chrm = chrm.strip()
    if not chrm.startswith("chr"):
        chrm = "chr"+chrm
    # coord should be in base 0 already, just need it as an int
    coord = int(coord)
    lmeth=int(lmeth)
    lunmeth=int(lunmeth)
    rmeth=int(rmeth)
    runmeth=int(runmeth)
    return (chrm,coord,lmeth,lunmeth,rmeth,runmeth)



    
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["combined=","leftErrorRate=","rightErrorRate="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    leftsuffix = None
    rightsuffix = None 
    
    for o, a in opts:
        if o=="--combined":
            infile = a
            print "combined", a
        elif o=="--leftErrorRate":
            leftErrorRate = float(a)
            print "Left Error Rate", a
        elif o=="--rightErrorRate":
            rightErrorRate = float(a)
            print "Right Error Rate", a

    
    assert infile != None
    
    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    
    binomOutput = csv.writer(open(infile+".binom","w"),delimiter='\t')
    

    for line in methfile:
        (chrm,coord,lmeth,lunmeth,rmeth,runmeth) = parseMethLine(line)
        
        binomOutput.writerow([chrm,coord,
                              lmeth,lunmeth,
                              rmeth,runmeth,
                              binom_test(lmeth,lmeth+lunmeth,leftErrorRate),
                              binom_test(rmeth,rmeth+runmeth,rightErrorRate)])