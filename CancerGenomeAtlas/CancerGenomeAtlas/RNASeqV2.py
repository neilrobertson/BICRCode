'''
Created on 29 Aug 2014

@author: neilrobertson
'''
from os import listdir
from os.path import isfile, join

from DataFile import DataFile
    
class RNASeqV2(object):
    def __init__(self, directory, outputFileName, rnaSeqFormat, cancerGenomeAtlas, delimiter = "\t"):
        from CancerGenomeAtlas import CancerGenomeAtlas
        self.directory = CancerGenomeAtlas.checkDirectory(directory)
        
        self.outputFileName = outputFileName
        self.rnaSeqFormat = rnaSeqFormat
        assert self.directory
        
        self.geneExpressionFiles = [f for f in listdir(self.directory) if (isfile(join(self.directory,f)) and f.endswith(rnaSeqFormat))]
        self.expressionDataMatrixDict = {}
        self.sampleIDs = []
        
        self.cancerGenomeAtlas = cancerGenomeAtlas
        self.outputDelimiter = delimiter
        
        print "A list of %s %s files has been found in directory %s" % (str(len(self.geneExpressionFiles)), self.rnaSeqFormat, self.directory)
        assert self.geneExpressionFiles
        
        
        
        
    def buildMatrix(self): 
        '''
        Builds a matrix of all the methylation files beta values using their genomic co-ords and generic data as a key
        '''
        from CancerGenomeAtlas import CancerGenomeAtlas
        
        if len(self.geneExpressionFiles) > 0:
            processedCounter = 1
            for fileName in self.geneExpressionFiles:
                fileName = CancerGenomeAtlas.checkFileName(fileName)
                if fileName:
                    fullFileName = self.directory + fileName
                    print "Processing file %s into the matrix!" % (str(processedCounter))   
                    
                    rnaFile = None
                    if self.rnaSeqFormat == r".rsem.genes.normalized_results":
                        rnaFile = RNASeqNormalizedGeneResultsFile(fullFileName)
                    elif self.rnaSeqFormat == r".rsem.genes.results":
                        rnaFile = RNASeqGeneResultsFile(fullFileName)   
                    
                    self.geneIDHeader, self.dataTypeHeader = rnaFile.openFileAndParseHeaders()
                    sampleCode = self.cancerGenomeAtlas.getSampleCodeFromFile(fileName)
                    if not sampleCode:
                        sampleCode = fileName
                    self.sampleIDs.append(sampleCode)
                    
                    for geneId, dataPoint in rnaFile.parseValues():
                        try:
                            dataPoints = self.expressionDataMatrixDict[geneId]
                            dataPoints.append(dataPoint)
                            self.expressionDataMatrixDict[geneId] = dataPoints
                        except:
                            self.expressionDataMatrixDict[geneId] = [dataPoint]
                    rnaFile.closeFile()
                    processedCounter += 1    
                    
                    
    def writeToOutput(self):
        '''
        Writes rnaSeq data matrix out to file
        '''
        outputFile = open(self.outputFileName,'w')
        delimitedSampleIDs = self.outputDelimiter.join(self.sampleIDs)
        outputFile.write("gene_id|gene_id_code" + self.outputDelimiter + delimitedSampleIDs + "\n")
        
        for key in self.expressionDataMatrixDict.keys():
            outputLine = ""
            expressionValues = self.expressionDataMatrixDict[key]
            
            outputLine = key + self.outputDelimiter
            for value in expressionValues:
                outputLine = outputLine + value.strip() + self.outputDelimiter
                
            if outputLine is not "":   
                outputLine = outputLine.strip(self.outputDelimiter) + "\n"
                outputFile.write(outputLine)
            
        outputFile.flush()
        outputFile.close()
        print "Successfully created rnaSeq matrix. Created file: %s" % (self.outputFileName)
        return self.outputFileName
    
    
    
    def dispose(self):
        self.expressionDataMatrixDict = None
           
            
    @staticmethod
    def getRNASequencers():
        return (r"UNC__IlluminaGA_RNASeqV2", r"UNC__IlluminaHiSeq_RNASeqV2")
    
    
    
    @staticmethod
    def getRNASeqFileSuffixes():
        return (r".junction_quantification.txt", r".rsem.genes.results", r".rsem.isoforms.results", 
                r".rsem.genes.normalized_results", r".rsem.isoforms.normalized_results", r".bt.exon_quantification.txt")
        
        
        
        
    @staticmethod
    def foldChange_RNASeqMatrix(fileName, annotationRows, arbitraryAddition = False, headers = True, delimiter = "\t", includeAllSiteAnnotation = False, logBase = 2):
        
        assert fileName
        assert int(annotationRows)
        print "Commencing RNA matrix log fold change filtering"
        outputFileName = ".".join(fileName.split(".")[:-1])+".geom_mean.foldChange_log%s.tsv" % str(logBase)
        
        matrixFile = open(fileName, 'r')
        output = open(outputFileName, 'w')
        
        for i, row in enumerate(matrixFile):
            if headers == True and i == 0:
                output.write(row)
            else:
                dataArray = row.strip().split(delimiter)[(annotationRows):]
                annotators = row.strip().split(delimiter)[:(annotationRows)]
                if annotators[annotationRows-1] in (delimiter, ""):
                    annotators[annotationRows-1] = "None"
                dataArray = RNASeqV2.foldChange_fromMean(dataArray, arbitraryAddition)
                failedRowCounter = 0
                if dataArray.any():
                    pyDataArray = dataArray.tolist()
                    if includeAllSiteAnnotation == False:
                        output.write(delimiter.join(annotators).strip(delimiter) + delimiter + delimiter.join(map(str, pyDataArray)).strip(delimiter) + "\n")
                else:
                    failedRowCounter += 1
        print "Totalled failed rows: %s" % (str(failedRowCounter))
        matrixFile.flush()
        matrixFile.close()
        output.flush()
        output.close()
        print "Completed RNA Seq filtering and log fold change from mean."
        return outputFileName



    @staticmethod
    def foldChange_fromMean(array, arbitraryAddition = False):
        import numpy as np
        import scipy.stats as scistat
        from Analysis.Filter import Filter
        array = Filter.convert_NAtoNaN(array)
        dataRow = np.array(array).astype(np.float)
        if Filter.test_allNaN(dataRow) == False:
            if Filter.test_containsNaN(dataRow) == True:
                dataRow = Filter.filterNonNumericCells_toZero(dataRow)
            if arbitraryAddition != False:
                assert float(arbitraryAddition)
                dataRow = RNASeqV2.addAbitraryAddition_nonZero(dataRow, arbitraryAddition)
            #mean = np.mean(dataRow, dtype=np.float64)
            mean = scistat.gmean(dataRow)
            for x in np.nditer(dataRow, op_flags=['readwrite']):
                x[...] = (np.log2(x) - np.log2(mean))
            return dataRow
        return None
    
    @staticmethod
    def foldChange_fromMedian(array, arbitraryAddition = False):
        import numpy as np
        from Analysis.Filter import Filter
        array = Filter.convert_NAtoNaN(array)
        dataRow = np.array(array).astype(np.float)
        if Filter.test_allNaN(dataRow) == False:
            if Filter.test_containsNaN(dataRow) == True:
                dataRow = Filter.filterNonNumericCells_toZero(dataRow)
            if arbitraryAddition != False:
                assert float(arbitraryAddition)
                dataRow = RNASeqV2.addAbitraryAddition_nonZero(dataRow, arbitraryAddition)
            median = np.median(dataRow)
            for x in np.nditer(dataRow, op_flags=['readwrite']):
                x[...] = (np.log2(x) - np.log2(median))
            return dataRow
        return None
    
    
    
    
    @staticmethod
    def addAbitraryAddition_nonZero(array, abitraryFigure):
        import numpy as np
        if type(array).__module__ == np.__name__:
            array = np.array(array).astype(np.float)
        assert float(abitraryFigure)
        for x in np.nditer(array, op_flags=['readwrite']):
            x[...] = x + abitraryFigure    
        return array
    
    @staticmethod
    def filterRNASeqData_standardDev(fileName, stdDev_cutoff = 2.0, headers = True, delimiter = "\t", includeAllSiteAnnotation = False):
        """
        Filters and cleans na values from data points and outputs data that has standard deviation above cutoff
        """
        from Analysis.Filter import Filter
        
        assert fileName
        assert float(stdDev_cutoff)
        print "Creating filtered methylation file. Cut-off Std Deviation: %s" % (str(stdDev_cutoff))
        
        outputFileName = ".".join(fileName.split(".")[:-1])+".filtered_stdDev%s.tsv" % str(stdDev_cutoff)
        
        matrixFile = open(fileName, 'r')
        output = open(outputFileName, 'w')
        
        minStdDev = 0
        maxStdDev = 0
        filteredCount = 0
        inclusionCount = 0
        
        for i, row in enumerate(matrixFile):
            if headers == True and i == 0:
                if includeAllSiteAnnotation == False:
                    sampleIDs = row.strip().split(delimiter)[1:]
                    headers = row.strip().split(delimiter)[0]
                    headerLine = headers.strip() + delimiter + delimiter.join(sampleIDs).strip(delimiter) + "\n"
                    output.write(headerLine)
                else:
                    output.write(row)
            else:
                dataArray = row.strip().split(delimiter)[1:]
                annotators = row.strip().split(delimiter)[0]
                dataArray, stdDev = Filter.filterMatrix_standardDev(dataArray)

                if stdDev:
                    if stdDev > maxStdDev:
                        maxStdDev = stdDev
                    if stdDev < minStdDev:
                        minStdDev = stdDev
                    if stdDev > stdDev_cutoff:
                        inclusionCount += 1
                        pyDataArray = dataArray.tolist()
                        if includeAllSiteAnnotation == False:
                            output.write(annotators.strip() + delimiter + delimiter.join(map(str, pyDataArray)).strip(delimiter) + "\n")
                        else:
                            output.write(annotators.strip(delimiter) + delimiter + delimiter.join(map(str, pyDataArray)).strip(delimiter) + "\n")
                    else:
                        filteredCount += 1
        print "Completed filtering RNASeq data. IncludedLines: %s Filtered Lines: %s" % (str(inclusionCount), str(filteredCount))
        print "Max StdDev: %s  Min StdDev: %s" % (maxStdDev, minStdDev)
        matrixFile.flush()
        matrixFile.close()
        output.flush()
        output.close()
        return outputFileName
    
    
    @staticmethod
    def rowScaleRawExpressionData(inputFilename, outputFilename, arbitraryAdditionn = float(0.001), containsHeader = True, annotationColumns = 1, delimiter = "\t"):
        '''
        Log^2 transforms raw RNA expression matrix and row scales each genes expression data into normalized scaled ranking between 0 and 1
        '''
        import numpy as np
        
        assert inputFilename
        assert outputFilename
        
        with open(inputFilename, "r") as inputFile:
            outputFile = open(outputFilename, "w")
            for i, row in enumerate(inputFile):
                if i == 0 and containsHeader == True:
                    outputFile.write(row)
                else:
                    rowParts = row.split(delimiter)
                    annotation = rowParts[annotationColumns - 1]
                    data = rowParts[annotationColumns:]
                    data = np.array(data).astype(np.float)           
                    for x in np.nditer(data, op_flags=['readwrite']):
                        x[...] = np.log2(x + arbitraryAdditionn)
                    absmin = abs(np.amin(data))
                    for x in np.nditer(data, op_flags=['readwrite']):
                        x[...] = x + absmin
                    absmax = abs(np.amax(data))
                    for x in np.nditer(data, op_flags=['readwrite']):
                        x[...] = x/absmax
                    pyDataArray = data.tolist()
                    outputFile.write(annotation.strip() + delimiter + delimiter.join(map(str, pyDataArray)).strip(delimiter) + "\n")
            outputFile.flush()
            outputFile.close()
        return outputFilename


    @staticmethod
    def filterGenes_byGeneIDs(inputFileName, outputFileName, containsHeader = True, delimiter = "\t"):
        '''
        Creates a filtered list of gene ids
        '''
        from DataFiles.DataFiles import DataFiles
        
        geneList = DataFiles.getGenesList()
        print geneList
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
                g_id = rowParts[0].strip().split("|")[0]
                print g_id
                if g_id in geneList:
                    identifiedGenes.append(row)
        
        if headerLine != None:
            outputFile.write(headerLine)
        for geneMutationData in identifiedGenes:
            outputFile.write(geneMutationData)
        print "Completed Gene ID list... Flushing resources."   
        inputFile.flush()
        inputFile.close()
        outputFile.flush()
        outputFile.close()
        return outputFileName


class RNASeqNormalizedGeneResultsFile(DataFile):
    def __init__(self, fileName, columnDelimiter = "\t"): 
        self.fileName = fileName
        self.rnaFile = None
        self.columnDelimiter = columnDelimiter
        assert self.fileName
        
        
        
        
    def openFileAndParseHeaders(self):
        ''' 
        Opens the RNA data file in TCGA format and returns header values
        '''
        self.rnaFile = open(self.fileName, 'r')
        assert self.rnaFile
        print "Successfully accessed file: %s" % (self.fileName)
            
        for i, line in enumerate(self.rnaFile):
            if i == 0:
                idHeader = line   
            elif i > 0:
                break
            
        idHeaderParts = idHeader.split(self.columnDelimiter)
        
        geneId = idHeaderParts[1]
        normalizedCount = idHeaderParts[0]
        return geneId, normalizedCount
    
    
    
    
    def parseValues(self):
        '''
        Parses all data lines in file and yields them as a generator
        '''
        for i, dataLine in enumerate(self.rnaFile):
            if i >= 1:
                lineParts = dataLine.split(self.columnDelimiter)
                
                geneID, normalizedCount = lineParts[0], lineParts[1]
                yield geneID, normalizedCount
            
             
             
                
    def closeFile(self):
        '''
        Closes methylation data file
        '''
        self.rnaFile.flush()
        self.rnaFile.close()   
        print "Successfully closed file: %s" % (self.fileName)
        print "=" *25      
    

    

class RNASeqGeneResultsFile(DataFile):
    def __init__(self, fileName, columnDelimiter = "\t"): 
        self.fileName = fileName
        self.rnaFile = None
        self.columnDelimiter = columnDelimiter
        assert self.fileName
        
        
        
    def openFileAndParseHeaders(self):
        ''' 
        Opens the RNA data file in TCGA format and returns header values
        '''
        self.rnaFile = open(self.fileName, 'r')
        assert self.rnaFile
        print "Successfully accessed file: %s" % (self.fileName)
            
        for i, line in enumerate(self.rnaFile):
            if i == 0:
                idHeader = line   
            elif i > 0:
                break
            
        idHeaderParts = idHeader.split(self.columnDelimiter)
        
        geneId = idHeaderParts[1]
        rawCount = idHeaderParts[0]
        
        scaledEstimate = idHeaderParts[1]  # @UnusedVariable
        transcriptId = idHeaderParts[2] # @UnusedVariable

        return geneId, rawCount
    
    
    def parseValues(self):
        '''
        Parses all data lines in file and yields them as a generator
        '''
        for i, dataLine in enumerate(self.rnaFile):
            if i >= 1:
                lineParts = dataLine.split(self.columnDelimiter)
                geneID, rawCount, scaledEstimate, transcriptID = lineParts[0], lineParts[1], lineParts[2], lineParts[3]  # @UnusedVariable
                yield geneID, rawCount
                
                
    def closeFile(self):
        '''
        Closes methylation data file
        '''
        self.rnaFile.flush()
        self.rnaFile.close()   
        print "Successfully closed file: %s" % (self.fileName)
        print "=" *25             



#RNASeqV2.filterGenes_byGeneIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/RNASeq/AllGenes.rsem.genes.normalized_results.output_matrix.foldChange_log2.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/CommonCRCMuts_vs_MethyltransferaseExpression/UNC__IlluminaGA_RNASeqV2.rsem.genes.results.output_matrix.Methyltransferases.tsv", containsHeader = True, delimiter = "\t")

#RNASeqV2.filterGenes_byGeneIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/RNASeq/AllGenes.rsem.genes.normalized_results.output_matrix.foldChange_log2.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNMT3B-RNA_BRAFv600E/UNC__IlluminaGA_RNASeqV2.rsem.genes.results.output_matrix.DNMT3B.tsv", containsHeader = True, delimiter = "\t")


#RNASeqV2.filterGenes_byGeneIDs("/mnt/50tb/privatedata/non-Adams/Brown-Borg.Holly/data/RNA_Seq/FPKM/FPKM_Matrix.csv", "/mnt/50tb/privatedata/non-Adams/Brown-Borg.Holly/data/RNA_Seq/FPKM/FPKM_Matrix_ChromatinGenes.csv", False)
#RNASeqV2.filterRNASeqData_standardDev("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/RNASeq_Clustering/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.normalized_results.output_matrix.foldChange_log2.intersectingIDs.SomMat_Meth450k.orderedClustered450kmeth.tsv",5)
#RNASeqV2.filterGenes_byGeneIDs("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/RNASeq_Clustering/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.normalized_results.output_matrix.foldChange_log2.intersectingIDs.SomMat_Meth450k.orderedClustered450kmeth.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/RNASeq_Clustering/Mann-Whitney/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.normalized_results.output_matrix.foldChange_log2.intersectingIDs.SomMat_Meth450k.orderedClustered450kmeth.MannWhitneyU_CIMP-HvsCIMP-L_fdr0.01.tsv")
#RNASeqV2.rowScaleRawExpressionData("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.normalized_results.output_matrix.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/RNAExpressionChanges/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.normalized_results.output_matrix.rowScaled.tsv")