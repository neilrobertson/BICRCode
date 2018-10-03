'''
Created on 29 Aug 2014

@author: neilrobertson
'''
import os, gc, time

from DNA_Methylation import DNA_Methylation
from RNASeqV2 import RNASeqV2
from Somatic_Mutations import Somatic_Mutations
from CNV_Low_Pass_DNASeq import CNV_Low_Pass_DNASeq
from Expression_Genes import Expression_Genes
from Expression_Protein import Expression_Protein


class CancerGenomeAtlas(object):
    def __init__(self, homeDirectory):
        self.homeDirectory = CancerGenomeAtlas.checkDirectory(homeDirectory)
        self.directoryContents = None
        self.outputPath = None
        
        self.fileManifestMap = None
        self.barcodeManifestMap = None
        
        self.fileDelimiter = "\t"
        self.availableDataTypes = {"DNA_Methylation":False, "Clinical":False, "RNASeqV2":False, "Somatic_Mutations":False, "CNV_SNP_Array":False, "Expression-Genes":False, "Expression-Protein":False}
        self.dataTypeOperations = {"Expression-Protein":self.buildProteinMatrix, "Expression-Genes":self.buildExpressionMatrix, "RNASeqV2":self.buildRNAMatrix, "DNA_Methylation":self.buildMethMatrix, "Somatic_Mutations":self.buildSomaticMutationsMatrix, "Clinical":self.performClinicalDataOperations, "CNV_Low_Pass_DNASeq":self.buildCNVSegmentationFiles}

        self.analyseFolderContents() 
        self.outputPath = self.createOutputFolder() 
        #self.buildRNAMatrix()
        self.performBuildingMatrixFromAvailableData()
        #self.buildCNVSegmentationFiles()
        print "[TCGA FILE PREPARATION COMPLETED!!]"
        
    
    def performClinicalDataOperations(self):
        pass
    
    def buildProteinMatrix(self):
        print "[CREATING PROTEIN EXPRESSION FILES]"
        for expressionFormat in Expression_Protein.getProteinExpressionFormatFolders():
            startTime = time.time()
            
            directory = self.homeDirectory + r"Expression-Protein/"+ expressionFormat + r"/Level_3"
            if CancerGenomeAtlas.checkDirectory(directory) != None:
                fileName = "Expression-Protein." + expressionFormat + ".output_matrix.tsv"

                outputDir = self.outputPath + r"/" + fileName
                
                proteinMatrix = Expression_Protein(directory, outputDir, self, expressionFormat)
                proteinMatrix.buildMatrix()
                buildTime = time.time()
                print "BuildTime: %s" % str(buildTime - startTime)
                
                outputFile = proteinMatrix.writeToOutput()
                endTime = time.time()
                proteinMatrix.dispose()
                print "EndTime %s" % str(endTime - startTime)
                return outputFile
        
        
    def buildExpressionMatrix(self):
        print "[CREATING GENE EXPRESSION FILES]"
        for expressionFormat in Expression_Genes.getExpressionFormatFolders():
            startTime = time.time()
            
            directory = self.homeDirectory + r"Expression-Genes/"+ expressionFormat + r"/Level_3"
            if CancerGenomeAtlas.checkDirectory(directory) != None:
                fileName = "Expression-Genes." + expressionFormat + ".output_matrix."
                if expressionFormat in ("UNC__AgilentG4502A_07_1", "UNC__AgilentG4502A_07_2", "UNC__AgilentG4502A_07_3"):
                    fileName += "lowessNormalised."
                fileName += "tsv"
                outputDir = self.outputPath + r"/" + fileName
                
                expressionMatrix = Expression_Genes(directory, outputDir, self, expressionFormat)
                expressionMatrix.buildMatrix()
                buildTime = time.time()
                print "BuildTime: %s" % str(buildTime - startTime)
                
                outputFile = expressionMatrix.writeToOutput()
                endTime = time.time()
                expressionMatrix.dispose()
                print "EndTime %s" % str(endTime - startTime)
                return outputFile
    
        
    def buildSomaticMutationsMatrix(self):
        print "[CREATING SOMATIC MUTATIONS MATRIX]"
        for mutationsSeqType in Somatic_Mutations.getSomaticMutationsSeqTypes():
            startTime = time.time()
            
            mutationsDir = self.homeDirectory + r"Somatic_Mutations/"+ mutationsSeqType + r"/Level_2"
            if CancerGenomeAtlas.checkDirectory(mutationsDir) != None:
                fileName = mutationsSeqType
                outputDir = self.outputPath + r"/" + fileName
                
                mutationsMatrix = Somatic_Mutations(mutationsDir, outputDir, mutationsSeqType, self)
                mutationsMatrix.run_MutSig()
                mutationsMatrix.buildMatrix()
                buildTime = time.time()
                print "BuildTime: %s" % str(buildTime - startTime)
                
                outputFile, mutationSignatureOutputFileName = mutationsMatrix.writeToOutput()
                endTime = time.time()
                mutationsMatrix.dispose()
                print "EndTime %s" % str(endTime - startTime) 
                ADDITIONAL_DATA_ANNOTATION = 10
                Somatic_Mutations.sortGenes_byMutationCount(mutationSignatureOutputFileName, ADDITIONAL_DATA_ANNOTATION)
                Somatic_Mutations.sortGenes_byMutationCount(outputFile)
                return outputFile
    
    
    def buildMethMatrix(self): 
        STD_DEV_FILTER_CUTOFF = 0.2
        
        print "[CREATING METHYFLATION FILES]"
        for methType in DNA_Methylation.getMethylationFormatFolders():
            startTime = time.time()
            
            methylationDir = self.homeDirectory + r"DNA_Methylation/"+ methType + r"/Level_3"
            if CancerGenomeAtlas.checkDirectory(methylationDir) != None:
                fileName = methType + ".output_matrix.tsv"
                outputDir = self.outputPath + r"/" + fileName
                
                methMatrix = DNA_Methylation(methylationDir, outputDir, self, methType)
                methMatrix.buildMatrix()
                buildTime = time.time()
                print "BuildTime: %s" % str(buildTime - startTime)
                
                outputFile = methMatrix.writeToOutput()
                endTime = time.time()
                methMatrix.dispose()
                print "EndTime %s" % str(endTime - startTime)
                DNA_Methylation.filterMethData_standardDev(outputFile, 5, STD_DEV_FILTER_CUTOFF)
                return outputFile
    
    
    def buildRNAMatrix(self):
        print "[CREATING RNA SEQ MATRIX FILES]"
        
        for rnaSeqType in RNASeqV2.getRNASequencers():
            startTime = time.time()
            
            rnaSeqDir = self.homeDirectory + r"RNASeqV2/"+ rnaSeqType + r"/Level_3"
            if CancerGenomeAtlas.checkDirectory(rnaSeqDir) != None:
                #fileTypes = [r".rsem.genes.results", r".rsem.genes.normalized_results"]
                fileTypes = [r".rsem.genes.normalized_results"]
                for fileType in fileTypes:
                    fileName = rnaSeqType + fileType + ".output_matrix.tsv"
                    print fileType
                    print "XX"

                    outputDir = self.outputPath + r"/" + fileName
                    
                    rnaSeqMatrix = RNASeqV2(rnaSeqDir, outputDir, fileType, self, delimiter = "\t")
                    rnaSeqMatrix.buildMatrix()
                    buildTime = time.time()
                    print "BuildTime: %s" % str(buildTime - startTime)
                    
                    outputFile = rnaSeqMatrix.writeToOutput()
                    endTime = time.time()
                    rnaSeqMatrix.dispose()
                    print "EndTime %s" % str(endTime - startTime)
                    if fileType == ".rsem.genes.normalized_results":
                        RNASeqV2.foldChange_RNASeqMatrix(outputFile, 1, 0.001)
                    return outputFile

        
    def buildCNVSegmentationFiles(self):
        #from ThirdPartyTools import GISTIC
        
        print "[CREATING CNV SEGMENTATION FILES AND ANALYZING]"
        for cnvFileType in CNV_Low_Pass_DNASeq.getCNVDataFolderTypes():
            startTime = time.time()
            
            cnvHomeDir = self.homeDirectory + r"CNV_SNP_Array/"+ cnvFileType + r"/Level_3"
            if CancerGenomeAtlas.checkDirectory(cnvHomeDir) != None:
                fileName = cnvFileType + ".globalSegmentationFile.tsv"

                outputFilename = self.outputPath + r"/" + fileName
                sampleListFilename = self.outputPath + r"/" + cnvFileType + r".sampleList.txt"
                print "Output Filename: %s" % outputFilename
                print "OutputSampleList Filename: %s" % sampleListFilename
                              
                cnvData = CNV_Low_Pass_DNASeq(cnvHomeDir, self.outputPath, self, delimiter = "\t")
                segmentationFilename, arraysListFile = cnvData.createSegmentationFile(outputFilename, sampleListFilename) #@UnusedVariable
                buildTime = time.time()
                print "BuildTime: %s" % str(buildTime - startTime)
                cnvData.dispose()
                return segmentationFilename
                #if segmentationFilename:
                #    gisticDir = self.createNewDirectory(self.outputPath + r"/GISTIC")
                #    if gisticDir != None:
                #        GISTIC.GISTIC.run_GISTIC(gisticDir, segmentationFilename, arraysListFile)
                #return segmentationFilename
            
        
    def performBuildingMatrixFromAvailableData(self):   
        for key in self.availableDataTypes.keys():
            if self.availableDataTypes[key] == True:
                self.dataTypeOperations[key]()
                gc.collect()
            
            
            
    def createOutputFolder(self):
        print "Checking for output directory"
        path = self.homeDirectory + r"Outputs"
        if not os.path.exists(path): 
            os.makedirs(path)
            print "Created new output directory"
        return CancerGenomeAtlas.checkDirectory(path)
        
        
    def createNewDirectory(self, directory):
        if not os.path.exists(directory): 
            try:
                os.makedirs(directory)
                return CancerGenomeAtlas.checkDirectory(directory)
            except:
                return None
        else:
            return directory
        
    def getOutputPath(self):
        return self.outputPath
    
    
    def analyseFolderContents(self):
        self.directoryContents = os.listdir(self.homeDirectory)
        for folder in self.directoryContents:
            if folder in self.availableDataTypes.keys():
                self.availableDataTypes[folder] = True
                print "Found data type in folder: %s" % (folder)
        print self.availableDataTypes 
        
        
        
    def getSampleCodeFromFile(self, fileName):
        if self.fileManifestMap == None:
            self.fileManifestMap = self.buildFileManifestDict()
        try:
            #return self.fileManifestMap[fileName]
            return "-".join(self.fileManifestMap[fileName].strip().split("-")[:4])[:-1]
        except:
            print "Not found sample id in fileName dictionary"
            return None
        
        
    def getSampleCodeFromBarcode(self, barcode):
        if self.barcodeManifestMap == None:
            self.barcodeManifestMap = self.buildBarcodeManifestDict()
        try:
            return self.barcodeManifestMap[barcode]
        except:
            print "Not found sample id in barcode dictionary. Barcode: %s" % (barcode)
            return None
      
        
    def buildBarcodeManifestDict(self):
        print "Building barcode manifest dictionary"
        manifestFileName = r"file_manifest.txt"
        manifestFileName = self.homeDirectory + manifestFileName
        
        manifestMap = {}
        
        manifestFile = open(manifestFileName, 'r')
        for i, line in enumerate(manifestFile):
            if i > 0:
                lineParts = line.split(self.fileDelimiter)
                platformType, center, platform, level, sample, barcode, fileName = lineParts[0], lineParts[1], lineParts[2], lineParts[3], lineParts[4], lineParts[5], lineParts[6]  # @UnusedVariable
                if sample != "selected_samples":
                    barcode = "-".join(barcode.strip().split("-")[:4])
                    manifestMap[barcode.strip()] = sample.strip()
        manifestFile.flush()
        manifestFile.close()
        return manifestMap
        
    
    
    def buildFileManifestDict(self):
        print "Building file manifest dictionary"
        #manifestFileName = r"file_manifest.txt"
        manifestFileName = r"FILE_SAMPLE_MAP.txt"
        manifestFileName = self.homeDirectory + manifestFileName
        
        manifestMap = {}
        
        manifestFile = open(manifestFileName, 'r')
        for i, line in enumerate(manifestFile):
            if i > 0:
                lineParts = line.strip().split(self.fileDelimiter)

                #platformType, center, platform, level, sampleId, barcode, fileName = lineParts[0], lineParts[1], lineParts[2], lineParts[3], lineParts[4], lineParts[5], lineParts[6]  # @UnusedVariable
                if len(lineParts) == 2:
                    fileName, sampleId = lineParts[0], lineParts[1]
                    if sampleId != "selected_samples":
                        manifestMap[fileName.strip()] = sampleId.strip()
        manifestFile.flush()
        manifestFile.close()
        return manifestMap
    
    
    
    @staticmethod
    def parseMethylationSampleIDs(sampleID):
        return "-".join(sampleID.split("-")[:3])
        #return sampleID
    @staticmethod
    
    def parseMethylationBarcode(sampleID):
        return "-".join(sampleID.split("-")[:4])
    
    @staticmethod
    def checkDirectory(directory):
        if os.path.isdir(directory):
            if directory[-1] != r"/":
                directory = directory + r"/"
            return directory
        else:
            return None
        
    @staticmethod
    def checkFileName(fileName):
        if fileName[0] == r".":
            return None
        return fileName
    
    
    

#tcga = CancerGenomeAtlas("/mnt/50tb/publicdata/TCGA/Skin_Cutaneous_Melanoma_28_08_2015/7d4e73d2-731f-4dc6-b873-bdaf0f788140/")
#cnv = CNV_Low_Pass_DNASeq("/mnt/50tb/publicdata/TCGA/Skin_Cutaneous_Melanoma_28_08_2015/7d4e73d2-731f-4dc6-b873-bdaf0f788140/CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3", "/mnt/50tb/privatedata/Rachael/TCGA_Melanoma/CNV/", tcga, delimiter = "\t")
#cnv.createSegmentationFile("/mnt/50tb/privatedata/Rachael/TCGA_Melanoma/CNV/Rachael.CRC.copyNumber.CNV_SNP_Array.SegmentationFile.tsv", "/mnt/50tb/privatedata/Rachael/TCGA_Melanoma/CNV/Rachael.CRC.copyNumber.CNV_SNP_Array.arraylist.txt")

#tcga = CancerGenomeAtlas("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_11_02_15/")
#cnv = CNV_Low_Pass_DNASeq("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_11_02_15/CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/CNV_Gistic/", tcga, delimiter = "\t")
#cnv.createSegmentationFile("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/CNV_Gistic/Douglas.CRC.copyNumber.BI__Genome_Wide.SegmentationFile.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/CNV_Gistic/Douglas.CRC.copyNumber.BI__Genome_Wide.arraylist.txt")

# mutationsMatrix = Somatic_Mutations("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_11_02_15/Somatic_Mutations/BCM__SOLiD_DNASeq/Level_2/", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_11_02_15/Outputs/BCM__SOLiD_DNASeq", "BCM__SOLiD_DNASeq", tcga)
# mutationsMatrix.buildMatrix()
# outputFile, mutationSignatureOutputFileName = mutationsMatrix.writeToOutput()
# mutationsMatrix.dispose()
# ADDITIONAL_DATA_ANNOTATION = 10
# Somatic_Mutations.sortGenes_byMutationCount(mutationSignatureOutputFileName, ADDITIONAL_DATA_ANNOTATION)
# Somatic_Mutations.sortGenes_byMutationCount(outputFile)