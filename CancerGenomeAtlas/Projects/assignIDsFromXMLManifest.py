from os import listdir
from os.path import isdir, join

#sudo pip install xmltodict
import xmltodict

manifestFilename = r"/mnt/50tb/privatedata/Rachael/TCGA_ECN_Restricted_Access/Output.xml"
dataRoot = r"/mnt/50tb/privatedata/Rachael/TCGA_ECN_Restricted_Access/data/FPKM"

manifestDict = {}
filenames = []

with open(manifestFilename, "r") as xmlFile:
    xmlManifest = xmltodict.parse(xmlFile.read())
    
    print "There are {0} results nodes in the xml manifest".format(str(len(xmlManifest['ResultSet']['Result'])))
    
    resultCounter = 0
    for resultNode in  xmlManifest['ResultSet']['Result']:
        filename = resultNode['files']['file']['filename']
        analysisID = resultNode['analysis_id']
        legacyID = resultNode['legacy_sample_id']
        
        filenames.append(filename)
        manifestDict[filename] = legacyID
        resultCounter += 1
        
assert len(filenames) == len(manifestDict.keys())
print "We have obtained {0} unique filenames from the xml manifest!".format(str(len(filenames)))

directoryContents = [x for x in listdir(dataRoot) if isdir(join(dataRoot,x))]
print "There are {0} data directories in the root.".format(str(len(directoryContents)))

processedFilenames = {}
for i in filenames:
    processed = i.strip().split(".")[2]
    processed = processed.strip().split("_")
    del processed[4]
    processed = "_".join(processed).strip()
    processedFilenames[processed] = i
    
matchedCount = 0
for dir in directoryContents:
    if dir in processedFilenames.keys():
        print "Matched xml manifest filename with directory!"
        sampleDir = dataRoot + "/" + dir + "/"
        processedSampleID = "-".join(manifestDict[processedFilenames[dir]].strip().split("-")[:4])[:-1]
        print processedSampleID
        with open(sampleDir+"sampleID.txt", "w") as sampleFile:
            sampleFile.write(processedSampleID.strip() + "\n")
        matchedCount += 1
        
print "Matched {0} directories with their sample ID!".format(str(matchedCount))