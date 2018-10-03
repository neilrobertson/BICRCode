'''
Created on 29 Aug 2014

@author: neilrobertson
'''

from os import listdir
from os.path import isfile, join

from DataFile import DataFile



class Expression_Genes(object):
    def __init__(self, directory, outputFileName, cancerGenomeAtlas, expressionFormat, delimiter = "\t"):
        from CancerGenomeAtlas import CancerGenomeAtlas
        
        self.directory = CancerGenomeAtlas.checkDirectory(directory)
        
        self.outputFileName = outputFileName
        self.expressionFormat = expressionFormat
        self.cancerGenomeAtlas = cancerGenomeAtlas
        assert self.directory
        
        self.geneExpressionFiles = [f for f in listdir(self.directory) if isfile(join(self.directory,f))]
        self.expressionDataMatrixDict = {}
        self.sampleIDs = []
        self.outputDelimiter = delimiter
        
        assert self.geneExpressionFiles
        print "A list of %s %s files has been found in directory %s" % (str(len(self.geneExpressionFiles)), self.expressionFormat, self.directory)


    def buildMatrix(self): 
        '''
        Builds a matrix of all the gene expression files using the gene name as a key
        '''
        from CancerGenomeAtlas import CancerGenomeAtlas
        
        if len(self.geneExpressionFiles) > 0:
            processedCounter = 1
            for fileName in self.geneExpressionFiles:
                fileName = CancerGenomeAtlas.checkFileName(fileName)
                if fileName:
                    fileName = self.directory + fileName
                    geneExpressionFile = GeneExpressionFile(fileName)
                    
                    sampleID, geneHeader, valueHeader = geneExpressionFile.openFileAndParseHeaders()  # @UnusedVariable
                    
                    sampleID = CancerGenomeAtlas.parseMethylationBarcode(sampleID)
                    print "Processing file %s into the matrix! SampleID: %s" % ((str(processedCounter)), sampleID)
                    if sampleID != None:
                        self.sampleIDs.append(sampleID)
                        
                        for  geneName, expressionValue in geneExpressionFile.parseValues():
                            
                            try:
                                expressionValues = self.expressionDataMatrixDict[geneName]
                                expressionValues.append(expressionValue)
                                self.expressionDataMatrixDict[geneName] = expressionValues
                            except:
                                self.expressionDataMatrixDict[geneName] = [expressionValue]
                    else:
                        print "No sampleID for file: %s" % (fileName)             
                    geneExpressionFile.closeFile()
                    processedCounter += 1
                   
                   
    def writeToOutput(self):
        '''
        Writes gene expression data matrix out to file
        '''
        outputFile = open(self.outputFileName,'w')

        delimitedSampleIDs = self.outputDelimiter.join(self.sampleIDs).strip()
        outputFile.write("Gene_Symbol" + self.outputDelimiter + delimitedSampleIDs + "\n")
        
        for geneName in self.expressionDataMatrixDict.keys():
            expressionValues = self.expressionDataMatrixDict[geneName]
            outputLine = geneName + self.outputDelimiter + self.outputDelimiter.join(expressionValues).strip() + "\n"
            outputFile.write(outputLine)
            
        outputFile.flush()
        outputFile.close()
        print "Successfully created gene expression matrix. Created file: %s" % (self.outputFileName)
        return self.outputFileName
    
    def dispose(self):
        self.expressionDataMatrixDict = None
    
    
    @staticmethod
    def getExpressionFormatFolders():
        return ("BI__HT_HG-U133A", "UNC__AgilentG4502A_07_1", "UNC__AgilentG4502A_07_2", "UNC__AgilentG4502A_07_3")     
        
        
class GeneExpressionFile(DataFile):
    def __init__(self, fileName, columnDelimiter = "\t"):
        self.fileName = fileName
        self.expressionFile = None
        self.columnDelimiter = columnDelimiter
        assert self.fileName
    
    def openFileAndParseHeaders(self):
        ''' 
        Opens the gene expression array file in TCGA format and returns header values
        '''
        self.expressionFile = open(self.fileName, 'r')
        assert self.expressionFile
        print "Successfully accessed file: %s" % (self.fileName)
        
        titleHeader = None
        idHeader = None
        for i, line in enumerate(self.expressionFile):
            if i == 0: idHeader = line
            elif i == 1: titleHeader = line    
            elif i > 1:
                break
            
        idHeaderParts = idHeader.split(self.columnDelimiter)
        titleHeaderParts = titleHeader.split(self.columnDelimiter)
        
        sampleID = idHeaderParts[1]
        geneHeader = titleHeaderParts[0]
        valueHeader = titleHeaderParts[1]
        
        return sampleID, geneHeader, valueHeader
    
    def parseValues(self):
        '''
        Parses all data lines in file and yields them as a generator
        '''
        for i, dataLine in enumerate(self.expressionFile):
            if i > 1:
                lineParts = dataLine.split(self.columnDelimiter)
                geneName, expressionValue = lineParts[0], lineParts[1]
                yield geneName, expressionValue
            
    def closeFile(self):
        '''
        Closes gene expression data file
        '''
        self.expressionFile.flush()
        self.expressionFile.close()   
        print "Successfully closed file: %s" % (self.fileName)
        print "=" *25 