'''
Created on 27 Nov 2014

@author: johncole
'''

import os
import getopt
import sys
from Data import ChIP as ChIPData
from Data import ChIPWithControl as ChIPWithControlData


#Method to deal with bed file entries:
def parseBedEntry (tokens):
    chrom = tokens[0]
    start = int(tokens[1])
    stop = int(tokens[2])
    assert start < stop, "End position of bed entry before start position, for line: " + tokens
    return chrom,start,stop


#Deals with the run options:
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["bedList=","fulloutPath=","summaryOutFile=","ChIP=","Control="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
        
    extend = None
    bedFiles = []
    fullOutPath = None
    summaryOutFile = None
    control = None
    ChIP = None
    dataSource = None
    usingControl = False
    
    for opt, o in opts:
        if opt == "--ChIP":
            ChIP = o
            print ChIP
        elif opt == "--Control":
            control = o
            usingControl = True
        elif opt == "--bedList":
            bedFiles.append(o)
        elif opt=="--fulloutPath":
            fullOutPath = o
        elif opt=="--summaryOutFile":
            summaryOutFile = open(o,'w')        
                
    assert len(bedFiles) > 0, "Need at least one --bedList"    
    assert ChIP != None, "No ChIP data"

print "### Loading Data ###"

#Loads the data:
if (usingControl):
    print "Using a control file"
    dataSource = ChIPWithControlData(ChIP,control,extend=extend)
    summaryOutFile.write("BedFile\tchipTags\tchipTagsRegionNormalised\tchipTagsLibraryNormalised\tchipTagsRegionLibraryNormalised\tcontrolTags\tcontrolTagsRegionNormalised\tcontrolTagsLibraryNormalised\tcontrolTagsRegionLibraryNormalised\tdiffTags\tdiffTagsRegionNormalised\tdiffTagsLibraryNormalised\tdiffTagsRegionLibraryNormalised\tratioTags\tratioTagsRegionNormalised\tratioTagsLibraryNormalised\tratioTagsRegionLibraryNormalised\tMeanChipTags\tMeanChipTagsRegionNormalised\tMeanChipTagsLibraryNormalised\tMeanChipTagsRegionLibraryNormalised\tMeanControlTags\tMeanControlTagsRegionNormalised\tMeanControlTagsLibraryNormalised\tMeanControlTagsRegionLibraryNormalised\tMeanDiffTags\tMeanDiffTagsRegionNormalised\tMeanDiffTagsLibraryNormalised\tMeanDiffTagsRegionLibraryNormalised\tMeanRatioTags\tMeanRatioTagsRegionNormalised\tMeanRatioTagsLibraryNormalised\tMeanRatioTagsRegionLibraryNormalised\n")
else:
    print "Not using a control file"
    dataSource = ChIPData(ChIP,extend=extend)
    summaryOutFile.write("BedFile\tTags\tTagsRegionNormalised\tTagsLibraryNormalised\tTagsRegionLibraryNormalised\tMeanTags\tMeanTagsRegionNormalised\tMeanTagsLibraryNormalised\tMeanTagsRegionLibraryNormalised\n")

print "### Processing Bed Files ###"

#Iterates through each bed file and gets the enrichment at each entry:
for bedFile in bedFiles:
    
    print "Processing file: " + str(fullOutPath + os.path.basename(bedFile))
    summaryTotals = []
    summaryCounts = []
    
    inFile = open(bedFile).readlines()
    fullOutFile = open(fullOutPath + os.path.basename(bedFile) + ".csv",'w') 
    
    #We print the headers:
    if (usingControl):
        fullOutFile.write("chrom\tstart\tstop\tchipTags\tchipTagsRegionNormalised\tchipTagsLibraryNormalised\tchipTagsRegionLibraryNormalised\tcontrolTags\tcontrolTagsRegionNormalised\tcontrolTagsLibraryNormalised\tcontrolTagsRegionLibraryNormalised\tdiffTags\tdiffTagsRegionNormalised\tdiffTagsLibraryNormalised\tdiffTagsRegionLibraryNormalised\tratioTags\tratioTagsRegionNormalised\tratioTagsLibraryNormalised\tratioTagsRegionLibraryNormalised\n")
    else:
        fullOutFile.write("chrom\tstart\tstop\tTags\tTagsRegionNormalised\tTagsLibraryNormalised\tTagsRegionLibraryNormalised\n")

    #We do the work:
    for entry in inFile:
        chrom,start,stop = parseBedEntry(entry.rstrip().split('\t'))
        
        if (usingControl):
            chipReadsDictionary,controlReadsDictionary = dataSource.getTagsInRegion(chrom, start, stop)
            normalisedTagsForEntry = dataSource.normaliseTags(chrom,start,stop,chipReadsDictionary,controlReadsDictionary)
        else:
            chipReadsDictionary = dataSource.getTagsInRegion(chrom, start, stop)
            normalisedTagsForEntry = dataSource.normaliseTags(chrom,start,stop,chipReadsDictionary)

        
        #Print the normalised tags for the entry:
        fullOutFile.write(str(chrom) + "\t" + str(start) + "\t" + str(stop))
        for i in normalisedTagsForEntry:
            fullOutFile.write("\t" + str(i))
        fullOutFile.write("\n")

        #Append the data to later calculate averages:
        if len(summaryCounts) == 0:
            for i in range(0,len(normalisedTagsForEntry)):
                summaryCounts.append(0)

        for i in range(0,len(summaryCounts)):
            if normalisedTagsForEntry[i] > 0.0:
                summaryCounts[i] = summaryCounts[i]+1
        
        if len(summaryTotals) == 0:
                summaryTotals = normalisedTagsForEntry
        else: 
            summaryTotals = [x + y for x, y in zip(summaryTotals,normalisedTagsForEntry)]
    
    
    #We print the total and average data for this bed file:
    summaryOutFile.write(os.path.basename(bedFile))    
    for i in summaryTotals:
        summaryOutFile.write("\t" + str(i))
    for i in range(0,len(summaryTotals)):
        if summaryCounts[i] > 0:
            summaryOutFile.write("\t" + str(float(summaryTotals[i]) / float(summaryCounts[i])))
        else:
            summaryOutFile.write("\t" + str(0.0))
    summaryOutFile.write("\n")
                             
    fullOutFile.close()

summaryOutFile.close()