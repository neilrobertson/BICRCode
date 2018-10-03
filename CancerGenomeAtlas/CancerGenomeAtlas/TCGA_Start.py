'''
Created on 16 Sep 2014

@author: neilrobertson

'''

import sys, getopt, time
from CancerGenomeAtlas import CancerGenomeAtlas

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["cancerGenomeAtlasDirectory="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
    cancerGenomeAtlasDirectory = None
    
    for o, a in opts:
        if o=="--cancerGenomeAtlasDirectory":
            cancerGenomeAtlasDirectory = a

    cancerGenomeAtlasDirectory = CancerGenomeAtlas.checkDirectory(cancerGenomeAtlasDirectory)            
    assert cancerGenomeAtlasDirectory != None
    
    if cancerGenomeAtlasDirectory[-1] != "/":
        cancerGenomeAtlasDirectory = cancerGenomeAtlasDirectory + "/"

    print "Initiating pre-processing of cancer genome atlas files..."
    startTime = time.time()
    CancerGenomeAtlas(cancerGenomeAtlasDirectory)
    
    print "Pre-processing completed! Time taken: %s" % (str(time.time() - startTime))
    
    