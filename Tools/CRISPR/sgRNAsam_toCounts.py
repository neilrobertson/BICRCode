'''
Created on 12 Aug 2015

@author: neilrobertson
'''

import getopt
import sys
import csv
from collections import defaultdict

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["sgRNALibrary=", "samFiles=", "sample-label=", "outputFilename="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
samFiles = None
sgRNAFilename = None 
sampleLabels = None
outputFilename = None
   
    
DELIMITER = ","
TAB_DEL = "\t"
    
for o, a in opts:
    if o=="--sgRNALibrary":
        sgRNAFilename = a
    elif o=="--samFiles":
        samFiles = a
    elif o=="--sample-label":
        sampleLabels = a
    elif o=="--outputFilename":
        outputFilename = a
        
assert samFiles != None
assert sgRNAFilename != None
assert sampleLabels != None
assert outputFilename != None

samFiles = [x.strip() for x in samFiles.split(",")]
sampleLabels = [x.strip() for x in sampleLabels.split(",")]

assert len(samFiles) == len(sampleLabels)

print "There are %s sam files provided for counting" % (str(len(samFiles)))

# Obtain dictionary of sgRNA names and related genes:
sgRNA_ids = []
sgRNA_id_dict = {}
with open(sgRNAFilename, "r") as sgRNAFile:
    for line in sgRNAFile:
        lineParts = line.strip().split(",")
        sgRNA_id = lineParts[0]
        gene_id = lineParts[2]
        sgRNA_ids.append(sgRNA_id)
        sgRNA_id_dict[sgRNA_id] = gene_id
        
print "We have a library of %s sgRNAs." % (str(len(sgRNA_ids)))

sgRNA_counts_dict = {}
for srRNA_id in sgRNA_ids:
    sgRNA_counts_dict[srRNA_id] = [0]*len(samFiles)
print "Initialised dictionary!"

fileCounter = 0
for samFilename in samFiles:
    
    print "*" * 25
    print "Performing counts on file: %s" % (samFilename)
    
    with open(samFilename, "r") as samFile:
        representedSequences = []
        testRNArepresentation = False
        
        unalignedCount, discoveredSeqs, totalSeqs = 0, 0, 0

        
        for line in samFile:
            lineParts = line.split(TAB_DEL)
            
            if lineParts[0] == "@HD":
                print "Read header, preparing to read seq lines..."
                
            elif lineParts[0] == "@SQ":
                representedSequences.append(lineParts[1].strip().split(":")[1])
            else:
                if not testRNArepresentation: 
                    nonRepresentedRNAs = [x for x in sgRNA_ids if x not in representedSequences]
                    totalRepresentedSeqs = len(representedSequences)
                    print "There are %s represented sequences. %s sequences not represented from sgRNA library." % (str(totalRepresentedSeqs), str(len(nonRepresentedRNAs)))
                    testRNArepresentation = True
                if lineParts[0] != "@PG":
                    if lineParts[2].strip() == "*":
                        unalignedCount += 1
                    else:
                        sgRNA_id = lineParts[2].strip()
                        sgRNA_counts_dict[sgRNA_id][fileCounter] += 1
                        discoveredSeqs += 1
                    totalSeqs += 1
        print "There was a %s percent sgRNA discovery rate for sgRNA sequences. With %s failures and %s found sgRNAs" % (str((float(discoveredSeqs)/float(totalSeqs))*100), str(unalignedCount), str(discoveredSeqs)) 
    fileCounter += 1
        

print "*"*25
print "Writing output...!"
with open(outputFilename, "w") as outputFile:
    
    output = csv.writer(outputFile, delimiter=TAB_DEL)
    
    header = ["sgRNA", "Gene"] + sampleLabels
    output.writerow(header)
    
    for sgRNA in sgRNA_ids:
        line = []
        line.append(sgRNA)
        line.append(sgRNA_id_dict[sgRNA])
        line = line + sgRNA_counts_dict[sgRNA]
        output.writerow(line)
print "Completed!"
        