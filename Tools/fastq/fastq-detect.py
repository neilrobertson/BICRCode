'''
Created on 10 Jul 2013

@author: mcbryan
'''

import getopt
import sys
from fastq.FastQ import FastQFile

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["fastq="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    fastq = None    
        
    for o, a in opts:
        if o=="--fastq":
            fastq = a
        
    assert fastq is not None
    
    firstread = None
    
    countmax = 10000
    
    minQual = None
    maxQual = None
    
    for count,read in enumerate(FastQFile(fastq)):
                
        if count == 0:
            firstread = read
        if count >= countmax:
            break
        
        minQual = read.minQual() if minQual is None else min(read.minQual(),minQual)
        maxQual = read.maxQual() if maxQual is None else max(read.maxQual(),maxQual)            
    
    # RANGES
    # Sanger -       Phred 33 (0   40)   : 33 - 73
    # Solexa -       Phred 64 (-5  40)   : 59 - 104
    # Illumina1.3 -  Phred 64 (0   40)   : 64 - 104
    # Illumina1.5 -  Phred 64 (0/2 40)   : 66 - 104 (0/1 scores not used)
    # Illumina1.8 -  Phred 33 (0   41)   : 33 - 74
    # Generic 64  -  Phred 64 (0   62)   : 64 - 126
    # Generic 33  -  Phred 33 (0   93)   : 33 - 126
    
    ranges = {"Sanger":(33,73),
              "Solexa":(59,104),
              "Illumina1.3":(64,104),
              "Illumina1.5":(65,104),
              "Illumina1.8":(33,74),
              "GenericPhred64":(64,126),
              "GenericPhred33":(33,126)}
    
    validRanges = set()
    for rangeid in ranges:
        low, high = ranges[rangeid]
        
        if minQual >= low and maxQual <= high:
            validRanges.add(rangeid)
            
    assert len(validRanges) > 0
    
    headerFormat = firstread.headerFormat()
    
    if headerFormat == "Unknown":
        remove = set()
    elif headerFormat == "Illumina1.8":
        remove = set(["Illumina1.3","Illumina1.5","Solexa","GenericPhred64"])
    elif headerFormat == "Illumina":
        remove = set(["Illumina1.8"])
    
    validRanges = validRanges - remove 
    
    assert len(validRanges) > 0
    
    # print most specific remaining in this order
    specificityOrder = ["Illumina1.5","Illumina1.3","Solexa","GenericPhred64","Sanger","Illumina1.8","GenericPhred33"]
    
    selected = None
    for specificity in specificityOrder:
        if specificity in validRanges:
            selected = specificity
            break
        
    assert selected != None
    
    outputMapping = {"Sanger":"phred33",
                     "Solexa":"solexa",
                     "Illumina1.3":"phred64",
                     "Illumina1.5":"phred64",
                     "Illumina1.8":"phred33",
                     "GenericPhred64":"phred64",
                     "GenericPhred33":"phred33"}
    
    print outputMapping[selected]