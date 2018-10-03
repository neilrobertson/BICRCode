"""

@author: Neil Robertson
"""

import sys, getopt, csv

DELIMITER = "\t"

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["inputFiles=", "consensusOutput="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

inputFiles = None
outputFilename = None

for o, a in opts:
    if o=="--inputFiles":
        inputFiles = a
    elif o=="--consensusOutput":
        outputFilename = a

assert inputFiles != None
assert outputFilename != None

inputFiles = inputFiles.strip().split(",")
print len(inputFiles)

positionCounterDict = {}
for filename in inputFiles:
    with open(filename, "r") as openedfile:
        bedfile = openedfile.readlines()
        for i, line in enumerate(bedfile):
            if i > 0:
                chrm, start = line.split(DELIMITER)[:2]
                try: positionCounterDict["%s:%s" % (chrm, start)] += 1
                except: positionCounterDict["%s:%s" % (chrm, start)] = 1
            
with open(outputFilename, "w") as outputFile:
    outputcsv = csv.writer(outputFile, delimiter=DELIMITER)
    for position in positionCounterDict.keys():
        chrm, start = position.split(":")
        count = positionCounterDict[position]
        outputcsv.writerow([chrm, start, count])