


class Analysis(object):
    def __init__(self):
        pass
    
    
    
    @staticmethod
    def rankSomaticMutations_byTreatmentContingency(inputRNAFilename, outputFilename, treatmentsCSVFilename, additionalHeaders = ("Entrez_Gene_Id","Chrom","Start_Position","End_Position","Strand","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2"), genesHeaderColumnName = "Hugo_Symbol", delimiter = "\t"):
        '''
        Contingency table approach for measuring significant difference in somatic mutations between different treatment groups.  
        The code builds a contingency table for each treatment pair using counts of mutated and non-mutated genes for each treatment.
        Use the additionalHeaders parameter to add additional headers that may occur in mutation signiture files.         
        '''
        from Utility.Utility import Utility
        import numpy as np
        import scipy.stats as scipystats
        
        USE_NUMPY_WRAPPING = True
        
        assert treatmentsCSVFilename
        assert inputRNAFilename
        assert outputFilename
        
        treatmentsDict = Analysis.buildTreatmentsLists(treatmentsCSVFilename, delimiter)

        if treatmentsDict not in ({},None):
            outputFile = open(outputFilename, "w")
            headerLine = genesHeaderColumnName + delimiter
            if additionalHeaders not in (None, ()):
                for addHeader in additionalHeaders:
                    headerLine += addHeader + delimiter
            for treat1 in treatmentsDict.keys():
                for treat2 in treatmentsDict.keys():
                    if treat1 != treat2:
                        headerLine += "%s_%s_pValue" % (treat1, treat2) + delimiter
            outputFile.write(headerLine.strip(delimiter) + "\n")
            
            with open(inputRNAFilename, "r") as inputRNAFile:
                columns, indexToName = Utility.getColumns(inputRNAFile) # @UnusedVariable

                geneList = columns[genesHeaderColumnName]
                
                totalMutations = len(geneList)
                counter = 0
                for i, gene in enumerate(geneList):
                    counter += 1
                    if counter % 1000 == 0:
                        print "Progress completed: %s percent.   Working on mutation at gene: %s" % (str((float(counter)/totalMutations)*100), gene)
                        
                    geneOutputLine = gene + delimiter
                    if additionalHeaders not in (None, ()):
                        for addHeader in additionalHeaders:
                            geneOutputLine += columns[addHeader][i] + delimiter
                    samplesDict = {}
                    for treatment in treatmentsDict.keys():
                        samplesDict[treatment] = []
                        samplesPerTreatment = treatmentsDict[treatment]
                        for sampleID in samplesPerTreatment:
                            try:
                                samplesDict[treatment].append(columns[sampleID][i])
                            except:
                                #print "Sample not found sampleID: %s" % (sampleID)
                                pass
                        if USE_NUMPY_WRAPPING:
                            samplesDict[treatment] = np.array(samplesDict[treatment]).astype(np.float)
                    for treat1 in treatmentsDict.keys():
                        for treat2 in treatmentsDict.keys():
                            if treat1 != treat2:
                                p = float(1.0)
                                try:
                                    p = Analysis.calculateContingencyP_expectedWeighted(samplesDict[treat1], samplesDict[treat2])
                                    geneOutputLine += str(p) + delimiter
                                except:
                                    geneOutputLine += str(p) + delimiter
                    outputFile.write(geneOutputLine.strip(delimiter) + "\n")
            outputFile.flush()
            outputFile.close()
            return outputFilename                    
    
    @staticmethod
    def calculateContingencyP_expectedWeighted(observed, expected):
        import scipy.stats as scipystats
        #Calculate the counts of mutation to non-mutant for each treatment
        observedMutationCount = Analysis.countNonZeros(observed)
        expectedMutationCount = Analysis.countNonZeros(expected)
        observedNonMutationCount = len(observed) - observedMutationCount
        expectedNonMutationCount = len(expected) - expectedMutationCount
        #Calculate percentile of the total in expected for mutant and non-mutant
        expectedMutationRatio = float(expectedMutationCount)/len(expected)
        expectedNonMutationRatio = float(expectedNonMutationCount)/len(expected)
        #Calculate the expected ratios of observed from the determined ratio in the expected figures
        predictedObservedMutation = len(observed)*expectedMutationRatio
        predictedObservedNonMutation = len(observed)*expectedNonMutationRatio
        try:
            #print [[observedMutationCount, observedNonMutationCount],[predictedObservedMutation, predictedObservedNonMutation]]
            oddsRatio, p = scipystats.fisher_exact([[observedMutationCount, observedNonMutationCount],[predictedObservedMutation, predictedObservedNonMutation]]) # @UnusedVariable
            return p
        except ValueError:
            return 1
    
    @staticmethod
    def countNonZeros(array):
        import numpy as np
        nonZeroCount = 0
        for x in np.nditer(array, op_flags=['readonly']):
            if x > 0: 
                nonZeroCount += 1
        return nonZeroCount
    
    @staticmethod
    def rankRNAExpressionTreatments_byMannWhitney(inputRNAFilename, outputFilename, treatmentsCSVFilename, genesHeaderColumnName = "gene_id|gene_id_code", delimiter = "\t"):
        '''
        Method uses Mann-Whitney test to determine variations in the distributions of RNA expression data between treatments.
        '''
        from Utility.Utility import Utility
        
        import numpy as np
        import scipy.stats as scistats
        
        TEST_TYPE = "mannwhitneyu"
        USE_NUMPY_WRAPPING = True
        
        assert treatmentsCSVFilename
        assert inputRNAFilename
        assert outputFilename
        
        treatmentsDict = Analysis.buildTreatmentsLists(treatmentsCSVFilename, delimiter)

        if treatmentsDict not in ({},None):
            outputFile = open(outputFilename, "w")
            headerLine = "gene_id" + delimiter
            for treat1 in treatmentsDict.keys():
                for treat2 in treatmentsDict.keys():
                    if treat1 != treat2:
                        headerLine += "%s_%s_U_%s" % (treat1, treat2, TEST_TYPE) + delimiter
                        headerLine += "%s_%s_p_%s" % (treat1, treat2, TEST_TYPE) + delimiter
            outputFile.write(headerLine.strip(delimiter) + "\n")
                        
            with open(inputRNAFilename, "r") as inputRNAFile:
                columns, indexToName = Utility.getColumns(inputRNAFile) # @UnusedVariable
                geneList = columns[genesHeaderColumnName]
                
                if geneList[0].find("|") != -1:
                    geneList = [x.split("|")[0] for x in geneList]
    
                for i, gene in enumerate(geneList):

                    geneOutputLine = gene + delimiter
                    samplesDict = {}
                    for treatment in treatmentsDict.keys():
                        samplesDict[treatment] = []
                        samplesPerTreatment = treatmentsDict[treatment]
                        for sampleID in samplesPerTreatment:
                            try:
                                samplesDict[treatment].append(columns[sampleID][i])
                            except:
                                #print "Sample not found sampleID: %s" % (sampleID)
                                pass
                        if USE_NUMPY_WRAPPING:
                            samplesDict[treatment] = np.array(samplesDict[treatment]).astype(np.float)
                            
                    for treat1 in treatmentsDict.keys():
                        for treat2 in treatmentsDict.keys():
                            if treat1 != treat2:
                                pValue = float(1.0)
                                U = float(0.0)
                                try:
                                    U, pValue = scistats.mannwhitneyu(samplesDict[treat1], samplesDict[treat2]) # @UnusedVariable
                                    geneOutputLine += str(U) + delimiter
                                    geneOutputLine += str(pValue) + delimiter
                                except:
                                    geneOutputLine += str(U) + delimiter
                                    geneOutputLine += str(pValue) + delimiter
                    outputFile.write(geneOutputLine.strip(delimiter) + "\n")
            
            outputFile.flush()
            outputFile.close()
            return outputFilename
        
         

    @staticmethod
    def buildTreatmentsLists(inputFilename, delimiter):
        from Utility.Utility import Utility
        assert inputFilename
        
        treatmentIDs = {}
        with open(inputFilename, "r") as inputFile:
            columns, indexToName = Utility.getColumns(inputFile)
            print "There are %s different treatment sets in the file: %s" % (str(len(indexToName.keys())), inputFilename)
            for index in indexToName.keys():
                columnHeader = indexToName[index]
                treatmentIDs[columnHeader] = tuple([x for x in columns[columnHeader] if x not in ("",None)])
                print "There are %s associated sample IDs in the treatment group labelled: %s" % (str(len(treatmentIDs[columnHeader])), columnHeader)
        return treatmentIDs
        
        
#Analysis.rankRNAExpressionTreatments_byMannWhitney("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/RNASeq_Clustering/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.normalized_results.output_matrix.foldChange_log2.intersectingIDs.SomMat_Meth450k.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/RNASeq_Clustering/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.normalized_results.output_matrix.foldChange_log2.intersectingIDs.SomMat_Meth450k.wilcoxTest_CIMPTreatments.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/CIMPClassificationSampleIDs.tsv")
#Analysis.rankSomaticMutations_byTreatmentContingency("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/BCM__IlluminaGA_DNASeq.mutationSignatureBySample.sorted.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Somatic_Mutations_Clustering/Somatic_Mutations_ContingencyComparisons/BCM__IlluminaGA_DNASeq.mutationSignatureBySample.ContingencyByTreatment.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/CIMPClassificationSampleIDs.tsv")