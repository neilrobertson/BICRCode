'''
Created on 22 May 2015

@author: neilrobertson
'''

import sys, getopt

DELIMITER = "\t"

def motif_block(datalist, n=2):
    for i in xrange(0, len(datalist), n):
        yield datalist[i:i+n]

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["motifInput=", "output=", "minSequenceLength=", "maxSequenceLength="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
        
inputFilename = None
outputFilename = None
minSequenceLength = 50
maxSequenceLength = 1000

for o, a in opts:
    if o=="--motifInput":
        inputFilename = a
    elif o=="--output":
        outputFilename = a
    elif o=="--minSequenceLength":
        minSequenceLength = int(a)
    elif o=="--maxSequenceLength":
        maxSequenceLength = int(a)

assert inputFilename != None
assert outputFilename != None

incorrectSeqLength = []
identifierDuplicates = []

motifDict = {}

print "Checking sequence lengths and unique identifiers..."
with open(inputFilename, "r") as inputFile:
    inputData = inputFile.readlines()
    for motifBlock in motif_block(inputData):
        identifier = motifBlock[0].strip()
        sequence = motifBlock[1].strip()
        if len(sequence) >= minSequenceLength and len(sequence) <= maxSequenceLength:
            try:
                priorSeq = motifDict[identifier]
                identifierDuplicates.append(identifier)
                if len(priorSeq) < len(sequence):
                    motifDict[identifier] = sequence
            except:
                motifDict[identifier] = sequence
        else:
            incorrectSeqLength.append(identifier)

print "Writing output..."
print "Sequence Count: %s" % str(len(motifDict.keys()))
with open(outputFilename, "w") as outputFile:
    for key in motifDict.keys():
        outputFile.write(key + "\n")
        outputFile.write(motifDict[key] + "\n")
        
print "Completed!!"
print "%s sequences rejected because of non-unique identifier." % (str(len(identifierDuplicates)))
print "%s sequences rejected due to length deficiencies." % (str(len(incorrectSeqLength)))