'''
Created on 15 Jan 2014

@author: johncole
'''
import getopt
import sys

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["input=","output=","mappings=","inputField=","mappingsLookupField=","mappingsOutputField="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
   
    tokens = []
    lookupDict = {}
    
    mappingsPrint = None
   
    for opt, a in opts:
        if opt=="--input":
            inFile = open(a).readlines()
        elif opt=="--output":
            outFile = open(a,'w')
        elif opt=="--mappings":
            mappings = open(a).readlines()
        elif opt=="--inputField":
            inputLookup = a
        elif opt=="--mappingsLookupField":
            mappingsLookup = a
        elif opt=="--mappingsOutputField":
            mappingsPrint = a 
   
    unmatched=0
    matched=0
   
    for line in mappings:
        
        tokens = line.rstrip().split('\t')
        #print tokens[int(mappingsLookup)].strip().lower()
        if mappingsPrint == None:
            #print tokens[int(mappingsLookup)]
            lookupDict[tokens[int(mappingsLookup)].strip().strip("\"")] = line.rstrip()
        else:
            lookupDict[tokens[int(mappingsLookup)].strip()] = tokens[int(mappingsPrint)].strip()
    #print lookupDict
    for line in inFile:
        tokens = line.rstrip().split('\t')
        #print tokens[int(inputLookup)]
	print tokens[int(inputLookup)].strip().replace("\"","")
        if lookupDict.has_key(tokens[int(inputLookup)].strip().replace("\"","")):
            
            outFile.write(line.rstrip() + "\t" + lookupDict[tokens[int(inputLookup)].strip().replace("\"","")] + "\n")
            matched += 1
            #outFile.write(lookupDict[tokens[int(inputLookup)].replace("\"","")] + "\n")

        else:
            #print tokens
            #print lookupDict.has_key(tokens[int(inputLookup)].strip().replace("\"",""))
            #print "Line did not match with mappings dictionary"
            unmatched += 1
            #outFile.write(line.rstrip() + "\n")
            #THIS HAS BEEN USED ALOT!!  ## outFile.write(line.rstrip() + "\t" + "NaN" + "\n")
            #outFile.write(line.rstrip() + "\tNA" + "\tNA" + "\tNA" + "\tNA" + "\tNA" + "\t\n")
       
    outFile.flush()
    outFile.close()
    
    print "Completed.  %s matched, %s unmatched." % (str(matched), str(unmatched))
