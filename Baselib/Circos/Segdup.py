

class SegdupCreator(object):
    def __init__(self):
        self.delimiter = "\t"

    def createFile(self, tool, inputFile, segdupFile, organismId, featureType, qualFilter, depthFilter, removeFeatures, intraInter, delimiter  = "\t"):
        
        if delimiter != self.delimiter: self.delimiter = delimiter
        
        lineProcessFactory = {"GASV":self.processLineGASV, "SVDetect":self.processLineSVDetect}
        
        ommittedLines = 0
        includedLines = 0
        for i, line in enumerate(inputFile):
            if i == 0: pass
            else:
                outputLine = lineProcessFactory[tool](line, organismId, featureType, qualFilter, depthFilter, removeFeatures, intraInter)
                if outputLine != None:
                    segdupFile.write(outputLine + "\n")
                    includedLines += 1
                else: ommittedLines += 1
        print "Filter has included %s lines and ommitted %s lines" % (str(includedLines), str(ommittedLines))
        return segdupFile
    
    
    
    def processLineGASV(self, line, organismId, featureType, qualFilter, depthFilter, removeFeatures, intraInter):
        lineParts=line.split(self.delimiter) #@UnusedVariable
        outputLine = None
        
        
    
        return outputLine
    
    
    
    def processLineSVDetect(self, line, organismId, featureType, qualFilter, depthFilter, removeFeatures, intraInter):
        lineParts=line.split(self.delimiter)
        
        chrType, svType, chr1, startEnd1, chr2, startEnd2, supportingPairs, totalScore, breakpoint1, breakpoint2 = lineParts[0], lineParts[1], lineParts[3], lineParts[4], lineParts[6], lineParts[7], int(lineParts[8]), float(lineParts[12]), lineParts[13], lineParts[14] #@UnusedVariable
        
        if intraInter != "ALL" and chrType.strip() != intraInter: return None
        
        if removeFeatures != None and svType in removeFeatures: return None
            
        if qualFilter != None and totalScore < qualFilter: return None
        
        if featureType != "ALL" and svType != featureType: return None

        if depthFilter != None and supportingPairs < depthFilter: return None

        if chr1.find("chr") != -1: chr1 = chr1.strip("chr")
        if chr2.find("chr") != -1: chr2 = chr2.strip("chr")
        
        chr1 = organismId + chr1
        chr2 = organismId + chr2
        
        return chr1 + " " + startEnd1.split("-")[0] + " " + startEnd1.split("-")[1] + " " + chr2 + " " + startEnd2.split("-")[0] + " " + startEnd2.split("-")[1]
    
    
    