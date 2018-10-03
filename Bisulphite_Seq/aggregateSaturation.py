import sys,getopt

def replicate_block(datalist, n=2):
    for i in xrange(0, len(datalist), n):
        yield datalist[i:i+n]

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["input=","output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

inputFilename=None
outputFilename=None

for o, a in opts:
    if o=="--input":
        inputFilename = a
    elif o=="--output":
        outputFilename = a

assert inputFilename != None
assert outputFilename != None

## Input file in bed format
with open(inputFilename, "r") as inputFile:
    
    with open(outputFilename, "w") as outputFile:
        
        for line in inputFile:
            chromPos = line.strip().split("\t")[:2]
            isRepresented = []
            for replicate in replicate_block(line.strip().split("\t")[2:]):
                if (int(replicate[0]) + int(replicate[1])) > 0:
                    isRepresented.append(1)
                else:
                    isRepresented.append(0)
            totalRepresentation = sum(isRepresented)
            
            outputFile.write(line.strip() + "\t" + str(totalRepresentation) + "\n")

        