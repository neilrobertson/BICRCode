import sys, getopt, os

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["inputDir=", "outputFilename=", "type="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
inputDir = None
outputFilename = None  
countsType = None 

    
DELIMITER = "\t"
    
for o, a in opts:
    if o=="--inputDir":
        inputDir = a.strip()
    elif o=="--outputFilename":
        outputFilename = a.strip()
    elif o=="--type":
        countsType = a.strip()
        
assert inputDir, outputFilename
assert countsType in ["htseq", "cuffdiff"]

print countsType


countDict = {}

if countsType == "htseq":
    
    print "Working on htseq counts directory files..."
    idList = []
    idsChecked = False
    files = next(os.walk(inputDir))[2]
    for thisfile in files:
        filename = thisfile.strip(".bam.HTSeq.txt")
        with open(inputDir+"/"+thisfile, "r") as currentFile:
            countsList = []
            print "Working on: {0}".format(thisfile)        
            for i, line in enumerate(currentFile):
                ensembl_id = line.strip().split()[0].strip()
                count = line.strip().split()[1].strip()
                
                countsList.append(count)
                if idsChecked:
                    if ensembl_id != idList[i]:
                        print("Failed....")
                        exit()
                else: 
                    idList.append(ensembl_id)
            countDict[filename] = countsList
    idsChecked = True
    
    print "Outputting matrix"

    with open(outputFilename, "w") as output:
        output.write("" + "\t" + "\t".join(idList) + "\n")
        for key in countDict.keys():
            output.write(key + "\t" + "\t".join(countDict[key]) + "\n")
    
    
elif countsType == "cuffdiff":
    
    print "Working on cuffdiff counts subdirectory files..."
    
    fpkm_dirs = next(os.walk(inputDir))[1]
    
    
    
    print fpkm_dirs
    for this_fpkm_dir in fpkm_dirs:
        filename = this_fpkm_dir.strip(".bam")
        idList = []
        with open(inputDir+"/"+this_fpkm_dir+"/"+"genes.fpkm_tracking", "r") as currentFile:
            localDict = {}
            print "Working on: {0}".format(this_fpkm_dir)
            for i, line in enumerate(currentFile):
                if i > 0:
                    ensembl_id = line.strip().split()[0].strip()
                    count = line.strip().split()[9].strip()
                    
                    localDict[ensembl_id] = count
                    idList.append(ensembl_id)
        countDict[filename] = localDict
        
    with open(outputFilename, "w") as output:
        
        #idList = set(idList).sort()
        print len(idList)
        idList = tuple(idList)
        outputHeader = "" + "\t" + "\t".join(idList) + "\n"
        output.write(outputHeader)
        
        for filename in countDict.keys():
            outputLine = [filename]
            for transcriptID in idList:
                outputLine.append(countDict[filename][transcriptID])
            output.write("\t".join(outputLine) + "\n")

print "Complete!"
            
        
                


    
