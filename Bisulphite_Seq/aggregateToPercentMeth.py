import csv,sys,getopt
import numpy as np

def replicate_block(datalist, n=2):
    for i in xrange(0, len(datalist), n):
        yield datalist[i:i+n]

delimiter = "\t"
threshold = 10
stdDevThreshold = None

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["file=","output=","threshold=","standardDevThreshold=","printStdDev="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
inputFilename = None
outputFilename = None
printStdDev = False

for o, a in opts:
    if o=="--file":
        inputFilename = a
    elif o=="--output":
        outputFilename = a
    elif o=="--threshold":
        threshold = int(a)
    elif o=="--standardDevThreshold":
        stdDevThreshold = float(a)
    elif o=="--printStdDev":
        printStdDev = True
        
assert inputFilename
assert outputFilename
        
outputFile = open(outputFilename,"w")
outputcsv = csv.writer(outputFile,delimiter=delimiter)
with open(inputFilename, "r") as inputFile:
    for line in inputFile:
        outputLine = []
        (chrm,position) = line.rstrip().split(delimiter)[:2]
        outputLine.append(chrm)
        outputLine.append(position)
        hasData = True # Test if all replicates have data before printing
        for replicateData in replicate_block(line.rstrip().split(delimiter)[2:]):
            meth, unmeth = float(replicateData[0]), float(replicateData[1])
            if (meth + unmeth) > threshold:
                outputLine.append(str((meth/(meth+unmeth))))
            else: hasData = False; break    
        if hasData == True: 
            if stdDevThreshold != None:
                stdDev = np.std(np.array(outputLine[2:]).astype(np.float))
                if printStdDev: outputLine.append(str(stdDev))
                if stdDev > stdDevThreshold: outputcsv.writerow(outputLine)
            else: outputcsv.writerow(outputLine)

outputFile.flush()
outputFile.close()