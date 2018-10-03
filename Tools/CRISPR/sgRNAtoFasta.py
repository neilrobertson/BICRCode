
import getopt
import sys

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["sgRNALibrary=", "fastaOut="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
fastaFilename = None
sgRNAFilename = None    
    
DELIMITER = ","
    
for o, a in opts:
    if o=="--sgRNALibrary":
        sgRNAFilename = a
    elif o=="--fastaOut":
        fastaFilename = a
        
numbers = [10,20,30,40,50,60]
import numpy as np
numbers = np.array(numbers).astype(np.float64)
average = np.average(numbers)
stdDev = np.std(numbers)
print average
print stdDev
for i in numbers:
    print (i - average)/stdDev

exit()

with open(fastaFilename, "w") as fastaFile:
    print "Writing to file: %s" % (fastaFilename)
    
    with open(sgRNAFilename, "r") as sgRNAFile:
        print "Loaded file: %s" % (sgRNAFilename) 
        
        for line in sgRNAFile:
            
            lineParts = line.strip().split(DELIMITER)
            sgRNA_id = lineParts[0].strip()
            sequence = lineParts[1].strip()
            geneName = lineParts[2].strip()
            
            fastaHeader = ">%s" % sgRNA_id
            
            fastaFile.write(fastaHeader + "\n")
            fastaFile.write(sequence + "\n")
            