"""
P.Adams lab.

Based on a concept from Landau et al; "Locally disordered methylation forms the basis for intra-tumor methylome 
variation in chronic lymphocytic leukemia"; Cancer Cell; Dec 2014

@author: Neil Robertson
"""

import getopt, sys, csv, gc, time, math

#Set file delimiter value
DELIMITER = "\t"

def parseContextLine(contextLine):
    contextLineParts = contextLine.split(DELIMITER)
    readID = contextLineParts[0]
    methStatus = contextLineParts[1]
    chrm = contextLineParts[2]
    position = contextLineParts[3]
    return readID, methStatus, chrm, int(position)


def parseRegionsLine(regionLine):
    regionChrm, regionStart, regionStop = regionLine.strip().split(DELIMITER)[:3]
    return regionChrm, int(regionStart), int(regionStop)

def addToPositionDictionary(position, readID, positionDict):
    try: positionDict[position].append(readID)
    except: positionDict[position] = [readID]
    return positionDict


def addToMethDictionary(methStatus, position, methDict):
    try: methDict[position].append(methStatus)
    except: methDict[position] = [methStatus]
    return methDict


def addToReadDictionary(methStatus, readID, readDict):
    try: readDict[readID].append(methStatus)
    except: readDict[readID] = [methStatus]
    return readDict


def calculateMethRatio(pos, methDict):
    calls = methDict[pos]
    meth = calls.count('+')
    unmeth = len(calls) - meth
    return meth, unmeth


def getProportionDiscordantReads(reads):
    discordantCount = 0
    for read in reads:
        if not all(read[0] == methCall for methCall in read):
            discordantCount += 1
    return discordantCount/float(len(reads)), discordantCount, (len(reads)-discordantCount)

def calculateGlobalPDM(readsDict):
    discordantCount = 0
    for r in readsDict.keys():
        if not all(readsDict[r][0] == methCall for methCall in readsDict[r]):
            discordantCount += 1
    return discordantCount/float(len(readsDict.keys()))


def getShannonEntropy(read):
    methCalls = read.count('+')
    methP = methCalls/float(len(read))
    unmethP = 1 - methP
    if methP == 1:
        shannonEnt = -(methP*math.log(methP, 2))
    elif unmethP == 1:
        shannonEnt = -(unmethP*math.log(unmethP, 2))
    else:
        shannonEnt = -((methP*math.log(methP, 2)) + (unmethP*math.log(unmethP, 2)))
    return shannonEnt


def getCumulativeShannonEntropy(reads):
    shannonEnts = []
    for r in reads:
        shannonEnts.append(getShannonEntropy(r))
    cumulativeEnt = sum(shannonEnts)/float(len(shannonEnts))
    return cumulativeEnt


def testReads(reads, readLengthThreshold):
    for r in reads:
        if len(r) < readLengthThreshold:
            reads.remove(r)
    if len(reads) >= 10:
        return True, reads
    else:
        return False, reads



def writePerSiteOutput(currentChrm, readDict, positionDict, methDict, outputcsv, readCoverageThreshold, readLengthThreshold):
    start = time.time()
    totalDiscordant = 0
    totalConcordant = 0
    orderedPositions = sorted(positionDict.keys())
    print "Creating output for chromosome: %s Positions: %s Total Reads: %s" % (currentChrm, str(len(orderedPositions)), str(len(readDict.keys())))
    for pos in orderedPositions:
        reads = positionDict[pos]
        if len(reads) >= readCoverageThreshold:
            readDetails = []
            for r in reads: readDetails.append(readDict[r])
            isRepresented, readDetails = testReads(readDetails, readLengthThreshold)
            if isRepresented:
                meth, unmeth = calculateMethRatio(pos, methDict)
                percentMeth = (meth/float(meth + unmeth))
                pdm, discordantCount, concordantCount = getProportionDiscordantReads(readDetails)
                totalDiscordant += discordantCount
                totalConcordant += concordantCount
                cumulativeShannonEnt = getCumulativeShannonEntropy(readDetails)
                outputLine = [currentChrm, str(pos), str(pos+1), str(percentMeth), str(concordantCount), str(discordantCount), str(len(reads)), str(pdm), str(cumulativeShannonEnt)]
                outputcsv.writerow(outputLine)
    print "Time taken: %s" % (time.time() - start)
    return totalDiscordant, totalConcordant


def writeWindowOutput(currentChrm, readDict, positionDict, methDict, outputcsv, windowSize, readCoverageThreshold, readLengthThreshold):
    start = time.time()
    totalDiscordant = 0
    totalConcordant = 0
    orderedPositions = sorted(positionDict.keys())
    startPosition = None
    for pos in orderedPositions:
        if startPosition == None:
            startPosition = pos
            
    print "Creating output for chromosome: %s Positions: %s Total Reads: %s" % (currentChrm, str(len(orderedPositions)), str(len(readDict.keys())))
    
    
    return totalDiscordant, totalConcordant


def writeRegionOutput(regionChrm, regionStart, regionStop, readDict, positionDict, methDict, outputcsv, readCoverageThreshold, readLengthThreshold):
    totalDiscordant, totalConcordant, totalMeth, totalUnmeth, representedPositions = [0,0,0,0,0]
    shannonEnts = []
    orderedPositions = sorted(positionDict.keys())
    totalPositions = len(orderedPositions)
    for pos in orderedPositions:
        reads = positionDict[pos]
        meth, unmeth = calculateMethRatio(pos, methDict)
        totalMeth += meth
        totalUnmeth += unmeth
        #if len(reads) >= readCoverageThreshold:
        readDetails = []
        for r in reads: readDetails.append(readDict[r])
        isRepresented, readDetails = testReads(readDetails, readLengthThreshold) #@UnusedVariable
        #if isRepresented:
        cumulativeShannonEnt = 0
        discordantCount, concordantCount = [0,0]
        if len(readDetails) > 0: 
            pdm, discordantCount, concordantCount = getProportionDiscordantReads(readDetails) #@UnusedVariable
            cumulativeShannonEnt = getCumulativeShannonEntropy(readDetails)
        totalDiscordant += discordantCount 
        totalConcordant += concordantCount
        if cumulativeShannonEnt:
            shannonEnts.append(cumulativeShannonEnt)
        representedPositions += 1
    regionPDM = 0
    regionShannonEnt = 0
    if totalDiscordant + totalConcordant > 0:
        percentMeth = (totalMeth/float(totalMeth + totalUnmeth))
        regionPDM = totalDiscordant/(float(totalDiscordant+totalConcordant))
        regionShannonEnt = sum(shannonEnts)/float(representedPositions)
        outputLine = [regionChrm, str(regionStart), str(regionStop), str(totalPositions), str(representedPositions), str(totalMeth), str(totalUnmeth), str(percentMeth), str(totalConcordant), str(totalDiscordant), str(regionPDM), str(regionShannonEnt)]
        outputcsv.writerow(outputLine)
    return None


def performPerSitePDM(contextFilename, outputFilename, readLengthThreshold, readCoverageThreshold):
    readDict = {}
    positionDict = {}
    methDict = {}
    totalConcordantReads = 0
    totalDiscordantReads = 0
    
    with open(outputFilename, "w") as outputFile:
        outputcsv = csv.writer(outputFile, delimiter=DELIMITER)
        
        outputcsv.writerow(["Chrm", "Start", "Stop", "PercentMeth", "Concordant", "Discordant", "OverlappingReads", "PDM", "ShannonEntropy"])
        
        currentChrm = None
        with open(contextFilename, "r") as contextFile:
            for i, line in enumerate(contextFile):
                readID, methStatus, chrm, position = (line)
                if currentChrm != chrm and currentChrm != None:
                    print "*" * 25
                    discordant, concordant = writePerSiteOutput(currentChrm, readDict, positionDict, methDict, outputcsv, readCoverageThreshold, readLengthThreshold)
                    totalConcordantReads += concordant
                    totalDiscordantReads += discordant
                    readDict = {}
                    positionDict = {}
                    methDict = {}         
                    gc.collect()
                    print "*" * 25
                currentChrm = chrm
                readDict = addToReadDictionary(methStatus, readID, readDict)
                positionDict = addToPositionDictionary(position, readID, positionDict)
                methDict = addToMethDictionary(methStatus, position, methDict)
                if (i % 1000000) == 0:
                    print "Currently working on line %s on chromosome %s." % (str(i), currentChrm)
    globalPDM = totalDiscordantReads / float(totalDiscordantReads + totalConcordantReads)
    print "Global PDM: %s Filename: %s" % (str(globalPDM), "/".join(contextFilename.split("/")[2:]))
    

def performWindowPDM(contextFilename, outputFilename, windowSize, readLengthThreshold, readCoverageThreshold):
    readDict = {}
    positionDict = {}
    methDict = {}
    totalConcordantReads = 0
    totalDiscordantReads = 0
    
    with open(outputFilename, "w") as outputFile:
        outputcsv = csv.writer(outputFile, delimiter=DELIMITER)
        
        outputcsv.writerow(["Chrm", "Start", "Stop", "PercentMeth", "Concordant", "Discordant", "OverlappingReads", "PDM", "ShannonEntropy"])
        
        currentChrm = None
        with open(contextFilename, "r") as contextFile:
            for i, line in enumerate(contextFile):
                readID, methStatus, chrm, position = parseContextLine(line)
                if currentChrm != chrm and currentChrm != None:
                    print "*" * 25
                    discordant, concordant = writeWindowOutput(currentChrm, readDict, positionDict, methDict, outputcsv, windowSize, readCoverageThreshold, readLengthThreshold)
                    totalConcordantReads += concordant
                    totalDiscordantReads += discordant
                    readDict = {}
                    positionDict = {}
                    methDict = {}         
                    gc.collect()
                    print "*" * 25
                currentChrm = chrm
                readDict = addToReadDictionary(methStatus, readID, readDict)
                positionDict = addToPositionDictionary(position, readID, positionDict)
                methDict = addToMethDictionary(methStatus, position, methDict)
                if (i % 1000000) == 0:
                    print "Currently working on line %s on chromosome %s." % (str(i), currentChrm)
    globalPDM = totalDiscordantReads / float(totalDiscordantReads + totalConcordantReads)
    print "Global PDM: %s Filename: %s" % (str(globalPDM), "/".join(contextFilename.split("/")[2:]))


def performRegionsPDM(contextFilename, outputFilename, regionsFilename, readLengthThreshold, readCoverageThreshold):
    
    with open(outputFilename, "w") as outputFile:
        outputcsv = csv.writer(outputFile, delimiter=DELIMITER)
        outputcsv.writerow(["Chrm", "Start", "Stop" "TotalCpGSites", "RepresentedSites", "Meth", "Unmeth", "PercentMeth", "TotalConcordant", "TotalDiscordant", "RegionPDM", "ShannonEntropy"])
    
        contextChrm = None
        contextPosition = None
        currentLine = None
        unfoundRegions = 0
        unfoundList = []
        
        with open(contextFilename, "r") as contextFile:
            regionsFile = open(regionsFilename, "r")
            regions = regionsFile.readlines()
            print "Working on %s regions..." % (len(regions))
            for i in range(0, len(regions)-1):
                regionChrm, regionStart, regionStop = parseRegionsLine(regions[i])
                print "Looking for region at %s : %s : %s" % (regionChrm, str(regionStart), str(regionStop))
                while contextChrm != regionChrm or contextPosition < regionStart:
                    currentLine = contextFile.readline()
                    readID, methStatus, contextChrm, contextPosition = parseContextLine(currentLine)
                if regionChrm == contextChrm and contextPosition >= regionStart and contextPosition <= regionStop:
                    print "Region found!"
                    readID, methStatus, chrm, position = parseContextLine(currentLine)   
                    readDict = {}
                    positionDict = {}
                    methDict = {}     
                    while position <= regionStop and regionChrm == chrm:
                        readDict = addToReadDictionary(methStatus, readID, readDict)
                        positionDict = addToPositionDictionary(position, readID, positionDict)
                        methDict = addToMethDictionary(methStatus, position, methDict)
                        currentLine = contextFile.readline()
                        readID, methStatus, chrm, position = parseContextLine(currentLine) 
                    writeRegionOutput(regionChrm, regionStart, regionStop, readDict, positionDict, methDict, outputcsv, readCoverageThreshold, readLengthThreshold)
                    print "Completed processing region."
                else:
                    print "Region not found"
                    unfoundRegions += 1
                    unfoundList.append("%s:%s" % (regionChrm, str(regionStart)))
            regionsFile.flush()
            regionsFile.close()
        print "Completed process... Undiscovered regions: %s" % (str(unfoundRegions))
        print unfoundList


#################################
######## BEGIN SCRIPT ###########
#################################

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["CpGContextFilename=","outputFilename=","regions=","windowSize=","readCoverageThreshold=","readLengthThreshold="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    
    contextFilename = None
    outputFilename =  None
    regionsFilename = None
    windowSize = None
    readLengthThreshold = 5
    readCoverageThreshold = 10
    
    for o, a in opts:
        if o=="--CpGContextFilename":
            contextFilename = a
        elif o=="--outputFilename":
            outputFilename = a
        elif o=="--regions":
            regionsFilename = a
        elif o=="--windowSize":
            windowSize = int(a)
        elif o=="--readCoverageThreshold":
            readCoverageThreshold = int(a)
        elif o=="--readLengthThreshold":
            readLengthThreshold = int(a)
            
    assert contextFilename != None
    assert outputFilename != None
    
    print "Input: %s" % (contextFilename)
    print "Output: %s" % (outputFilename)
            
    if regionsFilename:
        performRegionsPDM(contextFilename, outputFilename, regionsFilename, readLengthThreshold, readCoverageThreshold)
    elif windowSize:
        performWindowPDM(contextFilename, outputFilename, windowSize, readLengthThreshold, readCoverageThreshold)
    else:
        performPerSitePDM(contextFilename, outputFilename, readLengthThreshold, readCoverageThreshold)
    print "Process complete."    
