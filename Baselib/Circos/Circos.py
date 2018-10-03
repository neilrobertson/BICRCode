'''
Created on 1 Dec 2014

@author: neilrobertson
'''

'''
Created on 25 Jul 2014

@author: neilrobertson
'''

import getopt
import sys, time
from Segdup import SegdupCreator



if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["variance-file=","tool=","output-file=","qual-filter=","depth-filter=","organismId=","feature-type=","remove-features=", "intra-inter="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    inputFilename = None
    outputFilename= None
    tool = None
    qualFilter = None
    depthFilter = None
    featureType = None
    removeFeatures = None
    intraInter = None
    organismId = None
    
    FILE_DELIMITER = "\t"
    
    for o, a in opts:
        if o=="--variance-file":
            inputFilename = a
        elif o=="--tool":
            tool = a
        elif o=="--output-file":
            outputFilename = a
        elif o=="--qual-filter":    
            qualFilter = float(a)
        elif o=="--depth-filter":
            depthFilter = int(a)
        elif o=="--feature-type":
            featureType = a
        elif o=="--remove-features":
            removeFeatures = tuple(a.split(","))
        elif o=="--intra-inter":
            intraInter = a
        elif o=="--organismId":
            organismId = a
            
    assert inputFilename != None
    assert outputFilename != None
    assert tool.strip() in ("GASV", "SVDetect")
    assert organismId != None

    if featureType == None: featureType = "ALL"
    if intraInter == None: intraInter = "ALL"
    
    inputFile = open(inputFilename, "r")
    segdupFile = open(outputFilename, "w")
    
    startTime = time.time()
    
    print "Beginning to process and create segdup file..."
    segdupCreator = SegdupCreator()    
    segdupFile = segdupCreator.createFile(tool, inputFile, segdupFile, organismId, featureType, qualFilter, depthFilter, removeFeatures, intraInter, FILE_DELIMITER)
    
    endTime = time.time()
    segdupFile.flush()
    segdupFile.close()
    print "Completed creating Segdup file from format: %s" % (tool)
    print "Time taken: %s" % (str(endTime - startTime))