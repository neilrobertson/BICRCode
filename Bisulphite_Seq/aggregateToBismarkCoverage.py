'''
Created on 27 Feb 2015

@author: neilrobertson
'''

import getopt
import sys
import csv

def parseAggregateLine(line):
    (chrm,position,lmeth,lunmeth) = line.rstrip().split("\t")
    chrm = chrm.strip()
    if not chrm.startswith("chr"):
        chrm = "chr"+chrm
    # coord should be in base 0 already, just need it as an int
    position = int(position)
    meth=int(lmeth)
    unmeth=int(lunmeth)
    return (chrm,position,meth,unmeth)

def buildBismarkCovLine(chrm, position, meth, unmeth):
    percentMeth = 0.0
    if meth + unmeth > 0:
        percentMeth = (meth/(meth + float(unmeth)))*100 
    return [chrm, str(position+1), str(position+1), ("%.13f" % percentMeth), str(meth), str(unmeth)]

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["file=","output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
inputFilename = None
outputFilename = None

for o, a in opts:
    if o=="--file":
        inputFilename = a
    elif o=="--output":
        outputFilename = a
        
assert inputFilename != None
assert outputFilename != None

print "Input: %s" % (inputFilename)
print "Output: %s" % (outputFilename)
    
with open(inputFilename, "r") as inputFile:
    with open(outputFilename, "w") as outputFile:
        output = csv.writer(outputFile, delimiter="\t")
        for line in inputFile:
            (chrm, position, meth, unmeth) = parseAggregateLine(line)
            bismarkCovLine = buildBismarkCovLine(chrm, position, meth, unmeth)
            output.writerow(bismarkCovLine)
        