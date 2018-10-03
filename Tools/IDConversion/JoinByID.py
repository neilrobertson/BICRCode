'''
Created on 15 Jan 2014

@author: johncole
'''
import getopt
import sys

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["input=","output=","joinFeild=","printFeild="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    
    lookupDict = {}
    files = []
    joinFields = []
    printFeild = 0
    
    for opt, a in opts:
        if opt=="--input":
            files.append(a)
        elif opt=="--output":
            outFile = open(a,'w')
        elif opt=="--joinFeild":
            joinFields.append(int(a))
        elif opt=="--printFeild":
            printFeild = int(a)
            
    for path in files:
        inFile = open(path).readlines()
        print path
        
        for line in inFile:
            tokens = line.rstrip().split('\t')
            key = ""
            
            for i in joinFields:
                key += tokens[i] + "\t"
                
            
            if key not in lookupDict:
                temp = []
                lookupDict[key] = temp

            temp = lookupDict[key]
            temp.append(tokens[printFeild])
            lookupDict[key] = temp                
  
    for key in lookupDict:
        temp = lookupDict[key]
        outFile.write(key)
        for i in temp:
            outFile.write(i + "\t")
        outFile.write("\n")

    outFile.flush()
    outFile.close()