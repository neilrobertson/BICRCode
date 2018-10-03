'''
Created on 20 Aug 2014

@author: neilrobertson
'''


inputFile = 
outputFile = 
mappings = 
clockCSV = 

lookupDict = {}


for line in mappings:
    tokens = line.rstrip().split(',')
    lookupDict[tokens[int(0)]] = line.rstrip()


for line in inFile:
    tokens = line.rstrip().split('\t')
    cgID = tokens[0]
    if cgID.find("cg") != -1:
        if lookupDict.has_key(cgID):
            outFile.write(line.rstrip() + "\t" + lookupDict[tokens[int(inputLookup)].replace("\"","")] + "\n")
    else:
        outFile.write(line.rstrip() + "\t" + "\t" + "\n")
    
outFile.flush()
outFile.close()