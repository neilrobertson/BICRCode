'''
Created on 29 Aug 2014

@author: neilrobertson
'''

from os import listdir
from os.path import isfile, join

from DataFile import DataFile



class DNA_Methylation(object):
    def __init__(self, directory, outputFileName, cancerGenomeAtlas, dnaMethylationFormat, delimiter = "\t"):
        from CancerGenomeAtlas import CancerGenomeAtlas
        
        self.directory = CancerGenomeAtlas.checkDirectory(directory)
        
        self.outputFileName = outputFileName
        self.dnaMethylationFormat = dnaMethylationFormat
        self.cancerGenomeAtlas = cancerGenomeAtlas
        assert self.directory
        
        self.methylationFiles = [f for f in listdir(self.directory) if isfile(join(self.directory,f))]
        self.methylationDataMatrixDict = {}
        self.sampleIDs = []
        self.outputDelimiter = delimiter
        
        assert self.methylationFiles
        print "A list of %s %s files has been found in directory %s" % (str(len(self.methylationFiles)), self.dnaMethylationFormat, self.directory)
    
    
    
    def buildMatrix(self): 
        '''
        Builds a matrix of all the methylation files beta values using their genomic co-ords and generic data as a key
        '''
        from CancerGenomeAtlas import CancerGenomeAtlas
        
        if len(self.methylationFiles) > 0:
            processedCounter = 1
            for fileName in self.methylationFiles:
                fileName = CancerGenomeAtlas.checkFileName(fileName)
                if fileName:
                    fileName = self.directory + fileName
                    #while processedCounter < 6:
                    print "Processing file %s into the matrix!" % (str(processedCounter))    
                    methFile = MethylationArrayFile(fileName)
                    
                    barcode, cpgIDHeader, betaValueHeader, geneSymbolHeader, chromosomeHeader, genomicCoOrdHeader = methFile.openFileAndParseHeaders()  # @UnusedVariable
                    
                    sampleID = self.cancerGenomeAtlas.getSampleCodeFromBarcode(CancerGenomeAtlas.parseMethylationBarcode(barcode))
                    if sampleID != None:
                        self.sampleIDs.append(sampleID)
                        
                        for cpgId, betaValue, geneSymbol, chromosome, genomicCoOrd in methFile.parseValues():
                            key = ",".join([("chr"+str(chromosome)), genomicCoOrd, str(int(genomicCoOrd) + 1), cpgId, geneSymbol])
                            
                            try:
                                betaValues = self.methylationDataMatrixDict[key]
                                betaValues.append(betaValue)
                                self.methylationDataMatrixDict[key] = betaValues
                            except:
                                self.methylationDataMatrixDict[key] = [betaValue]
                    else:
                        print "No sampleID for barcode: %s" % (barcode)             
                    methFile.closeFile()
                    processedCounter += 1
      
             
             
                   
    def writeToOutput(self):
        '''
        Writes methylation data matrix out to file
        '''
        outputFile = open(self.outputFileName,'w')

        delimitedSampleIDs = self.outputDelimiter.join(self.sampleIDs)
        outputFile.write("Chromosome" + self.outputDelimiter + "start_Position" + self.outputDelimiter + "stop_Position" + self.outputDelimiter + "CpG_Ref" + self.outputDelimiter + "Gene_Symbol" + self.outputDelimiter + delimitedSampleIDs + "\n")
        
        for key in self.methylationDataMatrixDict.keys():
            outputLine = ""
            cpgKeyComponents = key.split(",")
            betaValues = self.methylationDataMatrixDict[key]
            
            for item in cpgKeyComponents:
                outputLine = outputLine + item.strip() + self.outputDelimiter
            for value in betaValues:
                outputLine = outputLine + value.strip() + self.outputDelimiter
                
            if outputLine is not "":   
                outputLine = outputLine.strip(self.outputDelimiter) + "\n"
                outputFile.write(outputLine)
            
        outputFile.flush()
        outputFile.close()
        print "Successfully created methylation matrix. Created file: %s" % (self.outputFileName)
        return self.outputFileName
    
    
    
    
    def dispose(self):
        self.methylationDataMatrixDict = None
        
   
   
    @staticmethod
    def getMethylationFormatFolders():
        return ("JHU_USC__HumanMethylation27","JHU_USC__HumanMethylation450")       
    
    
    
    @staticmethod
    def filterMethData_standardDev(fileName, annotationRows, stdDev_cutoff = 0.2, headers = True, delimiter = "\t", includeAllSiteAnnotation = False):
        """
        Filters and cleans na values from data points and outputs data that has standard deviation above cutoff
        """
        from Analysis.Filter import Filter
        
        assert fileName
        assert int(annotationRows)
        assert float(stdDev_cutoff)
        print "Creating filtered methylation file. Cut-off Std Deviation: %s" % (str(stdDev_cutoff))
        
        outputFileName = ".".join(fileName.split(".")[:-1])+".Filtered-StdDev%s.tsv" % str(stdDev_cutoff)
        
        matrixFile = open(fileName, 'r')
        output = open(outputFileName, 'w')
        
        for i, row in enumerate(matrixFile):
            if headers == True and i == 0:
                if includeAllSiteAnnotation == False:
                    sampleIDs = row.strip().split(delimiter)[(annotationRows):]
                    headers = row.strip().split(delimiter)[:(annotationRows)]
                    headerLine = headers[3].strip() + delimiter + delimiter.join(sampleIDs).strip(delimiter) + "\n"
                    output.write(headerLine)
                else:
                    output.write(row)
            else:
                dataArray = row.strip().split(delimiter)[(annotationRows):]
                annotators = row.strip().split(delimiter)[:(annotationRows)]
                if annotators[annotationRows-1] in (delimiter, ""):
                    annotators[annotationRows-1] = "None"
                dataArray, stdDev = Filter.filterMatrix_standardDev(dataArray)
                if stdDev and stdDev > stdDev_cutoff:
                    pyDataArray = dataArray.tolist()
                    if includeAllSiteAnnotation == False:
                        output.write(annotators[3].strip() + delimiter + delimiter.join(map(str, pyDataArray)).strip(delimiter) + "\n")
                    else:
                        output.write(delimiter.join(annotators).strip(delimiter) + delimiter + delimiter.join(map(str, pyDataArray)).strip(delimiter) + "\n")
        print "Completed filtering methylation data."
        matrixFile.flush()
        matrixFile.close()
        output.flush()
        output.close()
        return outputFileName
    
    @staticmethod
    def filterMethData_bySpecificSites(filename, filterAction, keyColumn = 1, delimiter = "\t", headers = True):
        """
        Filters methylation data based on CpG reference ids that are contained in the data files class.  The selection of 
        CpG references that will be filtered positively (maintained in file) or negatively (removed) is determined by the string: "filterAction"
        Current methods include:
        -SNPsRemoved
        -CpGIslandSites
        """

        assert filename
        assert int(keyColumn)
        
        filterSet, action, label = DNA_Methylation.getFilterMethods()[filterAction]()
        #print filterSet
        print label
        outputFileName = ".".join(filename.split(".")[:-1])+".%s.tsv" % str(label)
        outputFile = open(outputFileName, "w")
        
        sitesRemoved = 0
        with open(filename, 'r') as inputFile:
            for i, row in enumerate(inputFile):
                if headers == True and i == 0:
                    outputFile.write(row)
                else:
                    checkKey = row.split(delimiter)[keyColumn - 1].strip()
                    checkKey = checkKey.split("|")[0].strip() # Used in some instances where annotation coloum is concatenated
                    if (checkKey in filterSet) == action:
                        outputFile.write(row)
                    else:
                        sitesRemoved += 1
                        
        outputFile.flush()
        outputFile.close()
        print "Process of filtering %s has removed %s sites..." % (label, str(sitesRemoved))
        return outputFile
          
                    
    @staticmethod
    def getFilterMethods():
        """
        These filter methods correspond to the type of filtering that is required on the file which can be selected using the string keys
        """
        from DataFiles.DataFiles import DataFiles
        filterActionsDict = {"SNPsRemoved":DataFiles.getSNPSites_Meth450k, "CpGIslandSites":DataFiles.getCpGIslandSites_Meth450k}
        return filterActionsDict
   
    @staticmethod
    def matchCpGRefsToLocation(inputFilename, outputFilename, keyColumn, outputColumns, delimiter = "\t", header = True):
        """
        Re-pairs CpG reference ids to the location and gene data from the data file - CpGReferencesToMatch.txt
        """
        from DataFiles.DataFiles import DataFiles
        cpgRefs = DataFiles.getCpGReferences()
        cpgRefs = [x.split("|")[0].strip() for x in cpgRefs]
        print "Commencing mapping of CpG references to methylation location and gene data."
        outputFile = open(outputFilename, "w")
        
        with open(inputFilename) as inputFile:
            for i, row in enumerate(inputFile):
                if header == True and i == 0:
                    pass
                else:
                    rowParts = row.split(delimiter)
                    key = rowParts[keyColumn]
                    if key in cpgRefs:
                        outputLine = key + delimiter
                        for j in outputColumns:
                            outputLine += rowParts[j] + delimiter
                        outputFile.write(outputLine.strip() + "\n")
        outputFile.flush()
        outputFile.close()
        print "Completed mapping CpG references."
        return outputFilename

    @staticmethod
    def filterMethylationMatrix_byGeneIDs(inputFileName, outputFileName, containsHeader = True, delimiter = "\t"):
        '''
        Creates a filtered list of gene ids
        '''
        from DataFiles.DataFiles import DataFiles
        
        geneList = DataFiles.getGenesList()
        assert inputFileName
        assert len(geneList) > 0
        print "Beginning to filter list for corresponding gene IDs"
        inputFile = open(inputFileName, 'r')
        outputFile = open(outputFileName, 'w')
        
        headerLine = None
        identifiedGenes = []
        for i, row in enumerate(inputFile):
            if containsHeader == True and i == 0:
                headerLine = row
            else:
                rowParts = row.split(delimiter)
                g_id = rowParts[4].strip()   
                if g_id in geneList:
                    print "Appending site for gene: %s" % (g_id)
                    identifiedGenes.append(row)

        if headerLine != None:
            outputFile.write(headerLine)
        for geneMutationData in identifiedGenes:
            outputFile.write(geneMutationData)
        print "Completed Methylation matrix gene filtering... Flushing resources."   
        inputFile.flush()
        inputFile.close()
        outputFile.flush()
        outputFile.close()
        return outputFileName  
    
    @staticmethod
    def calculateMethylationStats(inputFilename, statsOutputFilename, sortedOutputFilename = None, annotationColumns = 1, headerLine = True, delimiter = "\t"):    
        import numpy as np
        from Utility.Utility import Utility
        from Analysis.Filter import Filter
        
        columns = {}
        columnHeaderDict = {}
        with open(inputFilename, "r") as inputFile:
            print "Opening file %s to parse per sample methylation data" % (inputFilename)
            columns, columnHeaderDict = Utility.getColumns(inputFile)
        
        ### CALCULATE STATS ###
        if len(columnHeaderDict) > annotationColumns:
            statsDict = {}
            meanDict = {}
            reverseDict = {}
            totalDict = {}
            
            for i in range(annotationColumns, len(columnHeaderDict.keys())):
                sampleID = columnHeaderDict[i]
                methArray = columns[sampleID]
                methArray = Filter.convert_NAtoNaN(methArray)
                dataPoints = np.array(methArray).astype(np.float)
                stdDev = np.std(dataPoints)
                mean = np.mean(dataPoints)
                median = np.median(dataPoints)
                total = np.sum(dataPoints)
                statsDict[sampleID] = (stdDev, mean, median, total)
                reverseDict[total] = sampleID  # This could be problematic if two sampleIDs share the same mean, although it seems relatively unlikely
                meanDict[sampleID] = mean
                totalDict[sampleID] = total
            
            
            outputSampleOrder = []
            if sortedOutputFilename:
                sortedTotals = sorted(totalDict.values(), reverse=True)
                for s_mean in sortedTotals:
                    outputSampleOrder.append(reverseDict[s_mean])
            else:
                for i in range(annotationColumns, len(columnHeaderDict.keys())):
                    outputSampleOrder.append(columnHeaderDict[i])
            
            usedStats = {"STD_DEV":0, "MEAN":1, "MEDIAN":2, "TOTAL":3}
            with open(statsOutputFilename, "w") as statsOutput:
                statsOutput.write("Statistic" + delimiter + delimiter.join(outputSampleOrder).strip() + "\n")
                for statType in usedStats.keys():
                    statLine = statType + delimiter
                    for s_id in outputSampleOrder:
                        statLine += str(statsDict[s_id][usedStats[statType]]) + delimiter
                    statsOutput.write(statLine.strip(delimiter) + "\n")
            
            if sortedOutputFilename:
                with open(sortedOutputFilename, "w") as outputFile:
                    outputFile.write("CpG-ID" + delimiter + delimiter.join(outputSampleOrder).strip() + "\n")
                    for j in range(0, len(columns[columnHeaderDict[0]])):
                        row = columns[columnHeaderDict[0]][j] + delimiter
                        for orderedID in outputSampleOrder:
                            row += columns[orderedID][j] + delimiter
                        outputFile.write(row.strip(delimiter) + "\n")
                print "Completed sorting and writing methylation matrix!"
            return sortedOutputFilename
                        
            

            
    @staticmethod
    def concatenateMatrixAnnotationColumns(inputFilename, outputFilename, removeNaN = True, delimiter = "\t", header = True):
        from Analysis.Filter import Filter
        import numpy as np
        
        with open(outputFilename, "w") as outputFile:
            print "Commencing annotation concatenation.  Will remove NA values: %s" % removeNaN
            inputFile = open(inputFilename, "r")
            for i, line in enumerate(inputFile):
                if header == True and i == 0:
                    sampleIDs = line.split(delimiter)[5:]
                    outputFile.write(line.split(delimiter)[3].strip() + " | " + line.split(delimiter)[4].strip() + delimiter + delimiter.join(sampleIDs).strip() + "\n")
                else:
                    isAllNaN = False
                    array = line.strip().split(delimiter)[5:]
                    dataPoints = Filter.convert_NAtoNaN(array)
                    dataPoints = np.array(dataPoints).astype(np.float)
                    isAllNaN = Filter.test_allNaN(dataPoints)
                    if isAllNaN == False:
                        pyDataArray = dataPoints.tolist()
                        outputFile.write(line.split(delimiter)[3].strip() + " | " + line.split(delimiter)[4].strip() + delimiter + delimiter.join(map(str, pyDataArray)).strip() + "\n")
            print "Completed. Flushing resources."
            inputFile.flush()
            inputFile.close()
        return outputFilename

    @staticmethod
    def removeNaN(inputFilename, outputFilename, removeNaN = True, delimiter = "\t", header = True):
        from Analysis.Filter import Filter
        import numpy as np
        
        with open(outputFilename, "w") as outputFile:
            print "Commencing annotation concatenation.  Will remove NA values: %s" % removeNaN
            inputFile = open(inputFilename, "r")
            for i, line in enumerate(inputFile):
                if header == True and i == 0:
                    outputFile.write(line.strip() + "\n")
                else:
                    isAllNaN = False
                    array = line.strip().split(delimiter)[5:]
                    dataPoints = Filter.convert_NAtoNaN(array)
                    dataPoints = np.array(dataPoints).astype(np.float)
                    isAllNaN = Filter.test_allNaN(dataPoints)
                    if isAllNaN == False:
                        pyDataArray = dataPoints.tolist()
                        outputFile.write(line.strip() + "\n")
            print "Completed. Flushing resources."
            inputFile.flush()
            inputFile.close()
        return outputFilename
            
        

class MethylationArrayFile(DataFile):
    def __init__(self, fileName, columnDelimiter = "\t"):
        self.fileName = fileName
        self.methylationFile = None
        self.columnDelimiter = columnDelimiter
        assert self.fileName
        
    
    
    
    def openFileAndParseHeaders(self):
        ''' 
        Opens the methylation array file in TCGA format and returns header values
        '''
        self.methylationFile = open(self.fileName, 'r')
        assert self.methylationFile
        print "Successfully accessed file: %s" % (self.fileName)
            
        for i, line in enumerate(self.methylationFile):
            if i == 0:
                idHeader = line
            elif i == 1:
                titleHeader = line    
            elif i > 1:
                break
            
        idHeaderParts = idHeader.split(self.columnDelimiter)
        titleHeaderParts = titleHeader.split(self.columnDelimiter)
        
        sampleID = idHeaderParts[1]
        cpgIDHeader = titleHeaderParts[0]
        betaValueHeader = titleHeaderParts[1]
        geneSymbolHeader = titleHeaderParts[2]
        chromosomeHeader = titleHeaderParts[3]
        genomicCoOrdHeader = titleHeaderParts[4]
        
        return sampleID, cpgIDHeader, betaValueHeader, geneSymbolHeader, chromosomeHeader, genomicCoOrdHeader
    
    
    
    
    def parseValues(self):
        '''
        Parses all data lines in file and yields them as a generator
        '''
        for i, dataLine in enumerate(self.methylationFile):
            if i > 1:
                lineParts = dataLine.split(self.columnDelimiter)
                cpgId, betaValue, geneSymbol, chromosome, genomicCoOrd = lineParts[0], lineParts[1], lineParts[2], lineParts[3], lineParts[4]
                yield cpgId, betaValue, geneSymbol, chromosome, genomicCoOrd
            
   
   
   
    def closeFile(self):
        '''
        Closes methylation data file
        '''
        self.methylationFile.flush()
        self.methylationFile.close()   
        print "Successfully closed file: %s" % (self.fileName)
        print "=" *25    


#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.CohortMatrices.CalorieRestriction.csv.RandomSubset", 1, 0.05, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.CohortMatrices.Exercise.csv.RandomSubset", 1, 0.05, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.CohortMatrices.WesternDiet.Old.csv.RandomSubset", 1, 0.05, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.CohortMatrices.WesternDiet.Young.csv.RandomSubset", 1, 0.05, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.ProcessedProbes.TCGA-CRC.NonTumour.DataMatrix.csv.RandomSubset", 1, 0.05, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)##

#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.CohortMatrices.CalorieRestriction.csv.RandomSubset", 1, 0.1, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.CohortMatrices.Exercise.csv.RandomSubset", 1, 0.1, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.CohortMatrices.WesternDiet.Old.csv.RandomSubset", 1, 0.1, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.CohortMatrices.WesternDiet.Young.csv.RandomSubset", 1, 0.1, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.ProcessedProbes.TCGA-CRC.NonTumour.DataMatrix.csv.RandomSubset", 1, 0.1, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)

#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.CohortMatrices.CalorieRestriction.csv.RandomSubset", 1, 0.15, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.CohortMatrices.Exercise.csv.RandomSubset", 1, 0.15, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.CohortMatrices.WesternDiet.Old.csv.RandomSubset", 1, 0.15, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.CohortMatrices.WesternDiet.Young.csv.RandomSubset", 1, 0.15, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/non-Adams/Luigi.Fontana/analysis/VariantProbes/cohortSize_4/Luigi-Fontana.CR.ProcessedProbes.TCGA-CRC.NonTumour.DataMatrix.csv.RandomSubset", 1, 0.15, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)


#DNA_Methylation.filterMethData_bySpecificSites("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/SecondPreProcessing/HumanMethylation450.NaN-Removed.tsv", "SNPsRemoved", keyColumn = 4, delimiter = "\t", headers = True)
#DNA_Methylation.removeNaN("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/PreProcessing/JHU_USC__HumanMethylation450.output_matrix.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/SecondPreProcessing/HumanMethylation450.output_matrix.NaNRemoved.tsv", removeNaN = True, delimiter = "\t", header = True)


#DNA_Methylation.calculateMethylationStats("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/PreProcessing/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.Filtered-StdDev0.2.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/GlobalMethylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.Filtered-StdDev0.2.Stats.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/GlobalMethylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.Filtered-StdDev0.2.MeanSorted.tsv", annotationColumns = 1, headerLine = True, delimiter = "\t")
#DNA_Methylation.calculateMethylationStats("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/PreProcessing/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/GlobalMethylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.Stats.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/GlobalMethylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.MeanSorted.tsv", annotationColumns = 1, headerLine = True, delimiter = "\t")
#DNA_Methylation.matchCpGRefsToLocation("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/PreProcessing/JHU_USC__HumanMethylation450.output_matrix.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/DNMT3B-CNV_Meth-Correlation/Douglas.CRC.copyNumber.DNMT3B.Meth450kData.SamplesOrdered.PearsonCorr0.4.bed", 3, [0,1,2,3,4], delimiter = "\t", header = True)
#DNA_Methylation.calculateMethylationStats("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.Filtered-StdDev0.2.CpGIslandSites.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.Filtered-StdDev0.2.CpGIslandSites.Stats.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.Filtered-StdDev0.2.CpGIslandSites.MeanSorted.tsv", annotationColumns = 1, headerLine = True, delimiter = "\t")
#DNA_Methylation.filterMethData_bySpecificSites("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.Filtered-StdDev0.2.tsv", "CpGIslandSites", keyColumn = 1, delimiter = "\t", headers = True)
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.tsv", int(1), stdDev_cutoff = 0.2, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.calculateMethylationStats("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/GlobalSampleStats/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.GlobalStats.tsv", sortedOutputFilename = None, annotationColumns = 1, headerLine = True, delimiter = "\t")
#DNA_Methylation.concatenateMatrixAnnotationColumns("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.output_matrix.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.NA-Removed.tsv")
#DNA_Methylation.filterMethData_bySpecificSites("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.NA-Removed.tsv", "SNPsRemoved", keyColumn = 1, delimiter = "\t", headers = True)

#DNA_Methylation.calculateMethylationStats("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_11_02_15/Outputs/JHU_USC__HumanMethylation450.output_matrix.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.output_matrix.BasicStats.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.SamplesMeanOrdered.tsv", 5)
#DNA_Methylation.filterMethData_bySpecificSites("/mnt/50tb/publicdata/TCGA/Prostatic_Adenocarcinoma_7_1_2015/Outputs/JHU_USC__HumanMethylation450.output_matrix.tsv", "CpGIslandSites",4)
#DNA_Methylation.calculateMethylationStats("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/JHU_USC__HumanMethylation450.output_matrix.filtered_lambda0.2.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Ordered_Methylation/JHU_USC__HumanMethylation450.output_matrix.filtered_lambda0.2.Statistics.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Ordered_Methylation/JHU_USC__HumanMethylation450.output_matrix.filtered_lambda0.2.SamplesOrdered.MEAN.tsv")
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Methylation_BoxPlots/JHU_USC__HumanMethylation450.output_matrix.IDsOnly.intersectingIDs.alignedToCluster.tsv", 1, 0.2, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
       
#DNA_Methylation.filterMethData_standardDev("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Douglas_LabMeetingImages/Meth450k_CIMPMarkerGenesOnly.intersectingIDs.alignedToCluster.CpGIslandSites.tsv", 1, 0.25, headers = True, delimiter = "\t", includeAllSiteAnnotation = True)
#DNA_Methylation.filterMethData_bySpecificSites("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Douglas_LabMeetingImages/Meth450k_CIMPMarkerGenesOnly.intersectingIDs.alignedToCluster.tsv", "CpGIslandSites")
#Filter on CpG Island sites        
#DNA_Methylation.filterMethylationMatrix_byGeneIDs("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/JHU_USC__HumanMethylation450.output_matrix.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Douglas_LabMeetingImages/Meth450k_CIMPMarkerGenesOnly.tsv")
#Add back cpg ref annotations
#DNA_Methylation.matchCpGRefsToLocation("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/JHU_USC__HumanMethylation450.output_matrix.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Methylation_Clustering/Meth450k_CpGCluster_CIMPRelatedCpGCluster.tsv", 3, (0,1,2,4))                