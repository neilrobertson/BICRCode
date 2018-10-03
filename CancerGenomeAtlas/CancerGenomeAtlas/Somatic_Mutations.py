'''
Created on 2 Sep 2014

@author: neilrobertson
'''
from os import listdir
from os.path import isfile, join
from collections import Counter

from DataFile import DataFile
from Utility.Utility import Utility

class Somatic_Mutations(object):
    def __init__(self, directory, outputFileName, fileFormat, cancerGenomeAtlas, delimiter = "\t"):
        from CancerGenomeAtlas import CancerGenomeAtlas
        print directory
        self.directory = CancerGenomeAtlas.checkDirectory(directory)
        
        self.outputFileName = outputFileName
        self.mutationsFormat = fileFormat
        self.cancerGenomeAtlas = cancerGenomeAtlas
        print self.directory
        assert self.directory
        
        self.mutationsFiles = [f for f in listdir(self.directory) if isfile(join(self.directory,f))]
        self.mutatedGenes = []
        self.fullMutationSignatures = []
        self.outputDelimiter = delimiter
        
        assert self.mutationsFiles
        print "A list of %s %s files has been found in directory %s" % (str(len(self.mutationsFiles)), self.mutationsFormat, self.directory)
    
    
    def run_MutSig(self):
        from ThirdPartyTools import MutSig
        from CancerGenomeAtlas import CancerGenomeAtlas
        
        if len(self.mutationsFiles) > 0:
            for fileName in self.mutationsFiles:
                fileName = CancerGenomeAtlas.checkFileName(fileName)
                if SomaticMutationsFile.isMafFile(fileName):
                    MutSig.MutSig.run_MutSig(self.directory + "/" + fileName, self.directory + "/" + ".".join(fileName.split(".")[:-1])+".MutSigResults")
         

    def buildMatrix(self):
        from CancerGenomeAtlas import CancerGenomeAtlas
        
        if len(self.mutationsFiles) > 0:
            processedCounter = 1
            for fileName in self.mutationsFiles:
                fileName = CancerGenomeAtlas.checkFileName(fileName)
                if SomaticMutationsFile.isMafFile(fileName):
                    fullFileName = self.directory + fileName
                    #while processedCounter < 6:
                    print "Processing file %s - %s" % (str(processedCounter), fileName)   
                    
                    somaticMutationsFile = SomaticMutationsFile(fullFileName, self.cancerGenomeAtlas)
                    self.uniqueIDs = somaticMutationsFile.getUniqueMutatedSampleIds()
                    self.uniqueSampleMutationsDict = {}
                    self.uniqueSampleFullMutationsDict = {}
                    for uniqueID in self.uniqueIDs:
                        self.uniqueSampleMutationsDict[uniqueID] = []
                        self.uniqueSampleFullMutationsDict[uniqueID] = []
                    print "File contains mutations for %s unique samples" % (str(len(self.uniqueIDs)))
                    
                    self.fileHeaders = tuple(somaticMutationsFile.openFileAndParseHeaders())

                    for geneID, entrezID, chrom, start, stop, strand, variantRef, refAllele, firstTumorSeqAllele, secondTumorSeqAllele, barcode in somaticMutationsFile.parseValues():
                        fullMutationSignature = self.outputDelimiter.join([geneID, entrezID, chrom, start, stop, strand, variantRef, refAllele, firstTumorSeqAllele, secondTumorSeqAllele]).strip(self.outputDelimiter)
                        barcode = "-".join(barcode.strip().split("-")[:4])
                        sampleID = self.cancerGenomeAtlas.getSampleCodeFromBarcode(barcode)
            
                        if geneID not in self.mutatedGenes:
                            self.mutatedGenes.append(geneID)
                        if fullMutationSignature not in self.fullMutationSignatures:
                            self.fullMutationSignatures.append(fullMutationSignature)
                        try:
                            self.uniqueSampleMutationsDict[sampleID].append(geneID)
                            self.uniqueSampleFullMutationsDict[sampleID].append(fullMutationSignature)
                        except:
                            print "HELP! My converted sample ids are not matching the correct tumor id... %s" % sampleID  
                            
                    somaticMutationsFile.closeFile()
                    processedCounter += 1 
    
    
    
    
    def writeToOutput(self):
        print "Writing paired mutations outputs file"
        geneOutputFileName = self.outputFileName + ".mutatedGeneBySample.tsv"
        geneOutputFile = open(geneOutputFileName,'w')

        geneOutputFile.write("Hugo_Symbol" + self.outputDelimiter + self.outputDelimiter.join(self.uniqueIDs).strip(self.outputDelimiter) + "\n")
        
        for mutatedGene in self.mutatedGenes:
            outputLine = ""
            outputLine += mutatedGene + self.outputDelimiter
            for sampleID in self.uniqueIDs:
                netMutationsList = self.uniqueSampleMutationsDict[sampleID]
                mutationCountAtGene = Counter(netMutationsList)[mutatedGene]
                outputLine += str(mutationCountAtGene) + self.outputDelimiter
                
            if outputLine is not "":   
                outputLine = outputLine.strip(self.outputDelimiter) + "\n"
                geneOutputFile.write(outputLine)
                
        geneOutputFile.flush()
        geneOutputFile.close()
        print "Successfully created somatic mutations matrix. Created file: %s" % (geneOutputFile) 
        
        geneSignatureOutputFileName = self.outputFileName + ".mutationSignatureBySample.tsv"
        geneSignatureOutputFile = open(geneSignatureOutputFileName,'w')
        geneSignatureOutputFile.write(self.outputDelimiter.join(self.fileHeaders).strip(self.outputDelimiter) + self.outputDelimiter + self.outputDelimiter.join(self.uniqueIDs).strip(self.outputDelimiter) + "\n")
        
        for mutatedGene in self.fullMutationSignatures:
            outputLine = ""
            outputLine += mutatedGene + self.outputDelimiter
            for sampleID in self.uniqueIDs:
                netMutationsList = self.uniqueSampleFullMutationsDict[sampleID]
                mutationCountAtGene = Counter(netMutationsList)[mutatedGene]
                outputLine += str(mutationCountAtGene) + self.outputDelimiter
                
            if outputLine is not "":   
                outputLine = outputLine.strip(self.outputDelimiter) + "\n"
                geneSignatureOutputFile.write(outputLine)
                
        geneSignatureOutputFile.flush()
        geneSignatureOutputFile.close()
        print "Successfully created somatic mutations matrix. Created file: %s" % (geneSignatureOutputFile) 
        return geneOutputFileName, geneSignatureOutputFileName
      
    
    def dispose(self):
        self.uniqueSampleFullMutationsDict = None
        self.fullMutationSignatures = None
        self.uniqueSampleMutationsDict = None
        self.mutatedGenes = None

    
    @staticmethod
    def getSomaticMutationsSeqTypes():
        return ("BCM__IlluminaGA_DNASeq", "BCM__SOLiD_DNASeq", "BCGSC__IlluminaHiSeq_DNASeq_automated", "BCM__IlluminaGA_DNASeq_automated", "BI__IlluminaGA_DNASeq_automated", "UCSC__IlluminaGA_DNASeq_automated", "BI__IlluminaGA_DNASeq_curated", "WUSM__IlluminaGA_DNASeq")
    

    @staticmethod
    def sortGenes_byMutationCount(inputFileName, annotationRows = 1, containsHeader = True, delimiter = "\t"):
        assert inputFileName
        print "Sorting somatic mutations file: %s" % inputFileName
        
        outFileName = ".".join(inputFileName.split(".")[:-1])+".sorted.tsv"
        
        inputFile = open(inputFileName, 'r')
        outputFile = open(outFileName, 'w')
        
        mutationsDict = {}
        headerLine = None
        for i, row in enumerate(inputFile):
            if containsHeader == True and i == 0:
                headerLine = row
            else:
                rowParts = row.split(delimiter)
                dataSet = rowParts[annotationRows:]    
                mutationTotal = 0
                for j in dataSet:
                    mutationTotal += int(j)
                mutationsDict[row] = mutationTotal       
                
        sortedL = sorted(mutationsDict.items(), key=lambda (k, v): v)
        
        if headerLine != None:
            outputFile.write(headerLine)
        for rowDict in reversed(sortedL):
            outputFile.write(rowDict[0])
        print "Completed sorting somatic mutations file.  Flushing resources..." 
        inputFile.flush()
        inputFile.close()
        outputFile.flush()
        outputFile.close()   
        
    
    @staticmethod
    def filterGenes_byGeneIDs(inputFileName, outputFileName, annotationColumns = 10, containsHeader = True, delimiter = "\t"):
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
                g_id = rowParts[0].strip()
                dataSet = rowParts[annotationColumns:]    
                if g_id in geneList:
                    identifiedGenes.append(row)
                    mutationTotal = 0
                    for j in dataSet:
                        mutationTotal += int(j)
                    print "Gene: %s Contains %s mutations." % (g_id, str(mutationTotal))
        
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
    
    
    @staticmethod
    def filterMafFile_bySampleIDs(inputFilename, outputFilename, delimiter = "\t", header = True, sampleColumnName = "Tumor_Sample_Barcode"):
        from DataFiles.DataFiles import DataFiles
        
        sampleIds = DataFiles.getSampleIDs_forFiltering()
        assert len(sampleIds) > 0
        #Ensure all sample Ids are the sample length, this is important later in method when they are compared to the .maf file.
        length = None
        for i in sampleIds:
            if length == None: length = len(i)
            if len(i) != length: assert "All sample IDs provided in \'SampleIDsFiltering.txt\' list not of identical length"
            
        assert inputFilename
        assert outputFilename
        
        print "Inputs correct.  GeneList compiled. Beginning to filter MAF file..."
        with open(inputFilename, "r") as inputFile:
            outputFile = open(outputFilename, "w")
            sampleColumnNumber = int(15)
            for i, line in enumerate(inputFile):
                if i == 0 and header == True:
                    outputFile.write(line)
                    #look up column where the correct sample IDs are located
                    headerParts = line.split(delimiter)
                    for j, name in enumerate(headerParts):
                        if name == sampleColumnName:
                            sampleColumnNumber = int(j)
                            print "Sample ID column found..."
                else:
                    currentSampleId = line.split(delimiter)[sampleColumnNumber]
                    if currentSampleId[:length] in sampleIds:
                        outputFile.write(line)

            outputFile.flush()
            outputFile.close()
        print "MAF File filtered successfully!"
        return outputFilename
                    
                    
                    
        
class SomaticMutationsFile(DataFile):
    def __init__(self, fileName, cancerGenomeAtlas, columnDelimiter = "\t"):
        self.fileName = fileName
        self.mutationsFile = None
        self.columnDelimiter = columnDelimiter
        self.cancerGenomeAtlas = cancerGenomeAtlas
        assert self.fileName
        
    
    def getUniqueMutatedSampleIds(self):
        columns, indexToName = self.buildColumns()  # @UnusedVariable
        ids = columns["Tumor_Sample_Barcode"]
        uniqueIDs = list()  
        map(lambda x: not x in uniqueIDs and uniqueIDs.append(x), ids)
        out = []
        failedCount = 0
        for i in uniqueIDs:
            barcode = "-".join(i.strip().split("-")[:4])
            sampleId = self.cancerGenomeAtlas.getSampleCodeFromBarcode(barcode)
            if sampleId != None:
                out.append(sampleId)
            else:
                failedCount += 1
        print "There were %s failed sampleID conversions from tumor sample barcode" % (str(failedCount))
        return tuple(out)
    
        
    def buildColumns(self):     
        columns, indexToName = Utility.getColumns(open(self.fileName, 'r'))
        return columns, indexToName
            
    
    def openFileAndParseHeaders(self):
        ''' 
        Opens the methylation array file in TCGA format and returns header values
        '''
        self.mutationsFile = open(self.fileName, 'r')
        assert self.mutationsFile
        print "Successfully accessed file: %s" % (self.mutationsFile)
        for i, line in enumerate(self.mutationsFile):
            if i == 0:
                idHeader = line    
            elif i > 0:
                break
        idHeaderParts = idHeader.split(self.columnDelimiter)
        hugoSymbol, entrezID, chrom, start, stop, strand, variantRef, refAllele, firstTumorSeqAllele, secondTumorSeqAllele  = Utility.unpack_fromPositions(idHeaderParts, (0,1,4,5,6,7,9,10,11,12))
        return hugoSymbol, entrezID, chrom, start, stop, strand, variantRef, refAllele, firstTumorSeqAllele, secondTumorSeqAllele

    
    def parseValues(self):
        '''
        Parses all data lines in file and yields them as a generator
        '''
        for i, dataLine in enumerate(self.mutationsFile):
            if i > 1:
                lineParts = dataLine.split(self.columnDelimiter)
                hugoSymbol, entrezID, chrom, start, stop, strand, variantRef, refAllele, firstTumorSeqAllele, secondTumorSeqAllele, barcode  = Utility.unpack_fromPositions(lineParts, (0,1,4,5,6,7,9,10,11,12,15))
                yield hugoSymbol, entrezID, chrom, start, stop, strand, variantRef, refAllele, firstTumorSeqAllele, secondTumorSeqAllele, barcode 
            
            
    def closeFile(self):
        '''
        Closes methylation data file
        '''
        self.mutationsFile.flush()
        self.mutationsFile.close()   
        print "Successfully closed file: %s" % (self.fileName)
        print "=" * 25  
        
    @staticmethod
    def isMafFile(filename):
        if filename.split(".")[-1] == "maf":
            return True
        else:
            return False    
        
#Somatic_Mutations.filterGenes_byGeneIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/SomaticMutations/AllSomaticMutations.mutatedGeneBySample.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/SomaticMutations/CRCCommonMutations.mutatedGeneBySample.tsv", annotationColumns = 1, containsHeader = True, delimiter = "\t")
#Somatic_Mutations.filterGenes_byGeneIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/SomaticMutations/AllSomaticMutations.mutationSignatureBySample.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNMT3B-Expression_vs_APCMutations/APCMutations.mutationSignatureBySample.tsv", annotationColumns = 10, containsHeader = True, delimiter = "\t")
#Somatic_Mutations.filterGenes_byGeneIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/SomaticMutations/AllSomaticMutations.mutationSignatureBySample.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNMT3B-Expression_vs_TP53Mutations/TP53Mutations.mutationSignatureBySample.tsv", annotationColumns = 10, containsHeader = True, delimiter = "\t")
#Somatic_Mutations.filterGenes_byGeneIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/SomaticMutations/AllSomaticMutations.mutatedGeneBySample.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/CommonCRCMuts_vs_MethyltransferaseExpression/BCM__IlluminaGA_DNASeq.mutationSignatureBySample.CommonCRCMutants.tsv", annotationColumns = 10, containsHeader = True, delimiter = "\t")
#Somatic_Mutations.filterGenes_byGeneIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/SomaticMutations/AllSomaticMutations.mutationSignatureBySample.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNMT3B-RNA_BRAFv600E/BCM__IlluminaGA_DNASeq.mutationSignatureBySample.BRAF.tsv", annotationColumns = 10, containsHeader = True, delimiter = "\t")


#Somatic_Mutations.filterGenes_byGeneIDs("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Somatic_Mutations_Clustering/BCM__IlluminaGA_DNASeq.mutatedGeneBySample.sorted.ordered450kMeth.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Douglas_LabMeetingImages/BCM__IlluminaGA_DNASeq.mutatedGeneBySample.sorted.ordered450kMeth.MethyltransferaseMutations.tsv")
#Somatic_Mutations.filterGenes_byGeneIDs("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/BCM__IlluminaGA_DNASeq.mutatedGeneBySample.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/MutSig/BCM__IlluminaGA_DNASeq.mutatedGeneBySample.MutSig-GloballySignificantlyMutatedGenes.tsv")        
#Somatic_Mutations.filterGenes_byGeneIDs("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Somatic_Mutations_Clustering/BCM__IlluminaGA_DNASeq.mutatedGeneBySample.sorted.ordered450kMeth.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Somatic_Mutations_Clustering/BCM__IlluminaGA_DNASeq.mutatedGeneBySample.sorted.ordered450kMeth.MutatableGenes.tsv")