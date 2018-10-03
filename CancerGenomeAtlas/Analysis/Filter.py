'''
Created on 3 Sep 2014

@author: neilrobertson
'''

import numpy as np

class Filter(object):
    def __init__(self):
        pass
    
    
    @staticmethod
    def filterfiles_forSharedSampleIDs(inputFilename, matchedIDFilename, outputFilename, delimiter = "\t"):
        from Utility.Utility import Utility
        
        matchedIDs = None
        with open(matchedIDFilename, 'r') as matchedIDsFile:
            for i, row in enumerate(matchedIDsFile):
                if i == 0: matchedIDs = tuple(row.strip().split(delimiter)[1:])
                else: break
                
        if matchedIDs:
            matchedIDs = tuple(list(set(matchedIDs)))
            print "Matched ID file found, containing %s IDs" % (str(len(matchedIDs)))
            outputFile = open(outputFilename, 'w')
            intersectingIDs = []

            with open(inputFilename) as inputFile:
                
                columns, indexToName = Utility.getColumns(inputFile)
                print "There are %s IDs in the input file" % (str(len(indexToName) - 1))

                intersectingIDs = [x for x in indexToName.values() if x in matchedIDs]
                intersectingIDs = tuple(list(set(intersectingIDs)))

                print "There have been %s intersecting ids discovered." % (str(len(intersectingIDs)))
                headerLine = indexToName[0].strip() + delimiter + indexToName[1].strip() + delimiter + indexToName[2].strip() + delimiter + indexToName[3].strip() + delimiter + indexToName[4].strip() + delimiter + delimiter.join(intersectingIDs).strip() + "\n"
                outputFile.write(headerLine)
                
                lineCount = 0
                for identifier in columns[indexToName[0]]:
                    line = identifier + delimiter + columns[indexToName[1]][lineCount]+ delimiter + columns[indexToName[2]][lineCount]+ delimiter + columns[indexToName[3]][lineCount]+ delimiter + columns[indexToName[4]][lineCount]+ delimiter
                    for consecutiveId in intersectingIDs:
                        line += str(columns[consecutiveId][lineCount]) + delimiter

                    outputFile.write(line.rstrip(delimiter) + "\n")
                    lineCount += 1
                outputFile.flush()
                outputFile.close()
                
    @staticmethod
    def getNonTumourSamples(inputFilename, outputFilename, sampleCode, delimiter = "\t", idDelimiter = "-"):
        from Utility.Utility import Utility
        with open(inputFilename) as inputFile:
            
            columns, indexToName = Utility.getColumns(inputFile)
            print "There are %s IDs in the input file" % (str(len(indexToName) - 1))
    
            selectedSamples = [x for x in indexToName.values() if x.split(idDelimiter)[-1] == sampleCode]
            print len(indexToName.values())
            print selectedSamples
            print "Total sample subset %s" % str(len(selectedSamples))
            
            
            with open(outputFilename, "w") as outputFile:

                headerLine = indexToName[0].strip() + delimiter + delimiter.join(selectedSamples).strip() + "\n"
                outputFile.write(headerLine)
                 
                lineCount = 0
                for identifier in columns[indexToName[0]]:
                    line = identifier + delimiter
                    for selectedSample in selectedSamples:
                        line += str(columns[selectedSample][lineCount]) + delimiter
         
                    outputFile.write(line.rstrip(delimiter) + "\n")
                    lineCount += 1

    
    
    
    @staticmethod
    def filterMatrix_standardDev(array):
        array = Filter.convert_NAtoNaN(array)
        #dataRow = np.array(rowParts[(annotationRows):]).astype(np.float)
        dataRow = np.array(array).astype(np.float)
        if Filter.test_allNaN(dataRow) == False:
            if Filter.test_containsNaN(dataRow) == True:
                dataRow = Filter.filterNonNumericCells_toMean(dataRow)
            stdDev = np.std(dataRow, dtype=np.float64)
            return dataRow, stdDev
        return None, None
                    
                    
    @staticmethod
    def convert_NAtoNaN(array):
        for i, point in enumerate(array):
            if point == "NA":
                array[i] = "NaN"
        return array
    
    
    @staticmethod         
    def filterNonNumericCells_toZero(row):
        if type(row).__module__ == np.__name__:
            row = np.array(row).astype(np.float)
        nanIndices = np.where(np.isnan(row))
        row[nanIndices] = np.float64(0.0)
        return row
    
               
    @staticmethod         
    def filterNonNumericCells_toMean(row):
        if type(row).__module__ == np.__name__:
            row = np.array(row).astype(np.float)
        nonNaNMean = row[~np.isnan(row)].mean()
        nanIndices = np.where(np.isnan(row))
        row[nanIndices] = nonNaNMean
        return row
    
    @staticmethod
    def test_containsNaN(array):
        if type(array).__module__ == np.__name__:
            array = np.array(array).astype(np.float)
        return np.isnan(array).any()
    
    @staticmethod
    def test_allNaN(array):
        if type(array).__module__ == np.__name__:
            array = np.array(array).astype(np.float)
        return np.isnan(array).all()
    
    @staticmethod
    def removeMatrixHeader(inputFilename, outputFilename, headerLines = 1, delimiter = "\t"):     
        with open(inputFilename, "r") as inputFile:
            outputFile = open(outputFilename, "w")
            for i, line in enumerate(inputFile):
                if i < headerLines: pass
                else:
                    outputFile.write(line)
            outputFile.flush()
            outputFile.close()
        return outputFilename


#Filter.getNonTumourSamples("/mnt/50tb/publicdata/TCGA/Breast_Cancer_14_10_15/Outputs/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.results.output_matrix.tsv", "/mnt/50tb/publicdata/TCGA/Breast_Cancer_14_10_15/Outputs/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.results.output_matrix.SolidTumourSamples.tsv", "01", delimiter = "\t", idDelimiter = "-")
#Filter.getNonTumourSamples("/mnt/50tb/publicdata/TCGA/Breast_Cancer_14_10_15/Outputs/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.results.output_matrix.tsv", "/mnt/50tb/publicdata/TCGA/Breast_Cancer_14_10_15/Outputs/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.results.output_matrix.SolidNonTumourSamples.tsv", "11", delimiter = "\t", idDelimiter = "-")
#Filter.getNonTumourSamples("/mnt/50tb/publicdata/TCGA/Breast_Cancer_14_10_15/Outputs/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.normalized_results.output_matrix.tsv", "/mnt/50tb/publicdata/TCGA/Breast_Cancer_14_10_15/Outputs/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.normalized_results.output_matrix.SolidTumourSamples.tsv", "01", delimiter = "\t", idDelimiter = "-")
#Filter.getNonTumourSamples("/mnt/50tb/publicdata/TCGA/Breast_Cancer_14_10_15/Outputs/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.normalized_results.output_matrix.tsv", "/mnt/50tb/publicdata/TCGA/Breast_Cancer_14_10_15/Outputs/UNC__IlluminaHiSeq_RNASeqV2.rsem.genes.normalized_results.output_matrix.SolidNonTumourSamples.tsv", "11", delimiter = "\t", idDelimiter = "-")

#Filter.filterfiles_forSharedSampleIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/SecondPreProcessing/HumanMethylation450.NaN-Removed.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNMT3B-RNA_BRAFv600E/BRAFMutations.DNMT3BExpression.MergedMatrix.csv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/SecondPreProcessing/HumanMethylation450.NaN-Removed.SharedSamples.tsv", delimiter = "\t")
#Filter.filterfiles_forSharedSampleIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/CNV_Gistic/Douglas.CRC.copyNumber.BI__Genome_Wide.SegmentationFile.sorted.DNMT3B.sorted.unique.cnvWithID.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/PreProcessing/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.Filtered-StdDev0.2.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/DNMT3B-CNV_Meth-Correlation/Douglas.CRC.copyNumber.DNMT3B.Meth450kSharedSamples.tsv", delimiter = "\t")
#Filter.filterfiles_forSharedSampleIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/SomaticMutations/CRCCommonMutations.mutatedGeneBySample.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/Clustering/MutationsHeatmaps/HumanMethylation450k.NoNoise.TumourSample.ClusteredSampleOrder.csv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/Clustering/MutationsHeatmaps/CRCCommonMutations.mutatedGeneBySample.MethClusterSharedSamples.tsv", delimiter = "\t")

#Filter.getNonTumourSamples("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/RNASeq/AllGenes.rsem.genes.results.output_matrix.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/RNASeq/AllGenes.rsem.genes.results.output_matrix.TumourSamples.tsv", "01", delimiter = "\t", idDelimiter = "-")
#Filter.getNonTumourSamples("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/RNASeq/AllGenes.rsem.genes.results.output_matrix.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/RNASeq/AllGenes.rsem.genes.results.output_matrix.NonTumourSamples.tsv", "11", delimiter = "\t", idDelimiter = "-")
#Filter.filterfiles_forSharedSampleIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNMT3B-RNA_BRAFv600E/BRAFMutations.DNMT3BExpression.MergedMatrix.csv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.tsv", delimiter = "\t")
#Filter.getNonTumourSamples("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.tsv", "01", delimiter = "\t", idDelimiter = ".")
#Filter.filterfiles_forSharedSampleIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.output_matrix.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNMT3B-RNA_BRAFv600E/BRAFMutations.DNMT3BExpression.MergedMatrix.csv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/JHU_USC__HumanMethylation450.SharedSamples.tsv", delimiter = "\t") 
    
#Filter.removeMatrixHeader("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Methylation_BoxPlots/JHU_USC__HumanMethylation450.output_matrix.IDsOnly.intersectingIDs.alignedToCluster.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Methylation_BoxPlots/JHU_USC__HumanMethylation450.output_matrix.IDsOnly.intersectingIDs.alignedToCluster.noheader.tsv")       
#Filter.filterfiles_forSharedSampleIDs("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/JHU_USC__HumanMethylation450.output_matrix.IDsOnly.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Methylation_Clustering/JHU_USC__HumanMethylation450.output_matrix.filtered_lambda0.2.CpGIslandSites.intersectingIDs.SomMat_IlluminaGA.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Methylation_BoxPlots/JHU_USC__HumanMethylation450.output_matrix.IDsOnly.intersectingIDs.tsv") 
#Filter.filterfiles_forSharedSampleIDs("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Douglas_LabMeetingImages/Meth450k_CIMPMarkerGenesOnly.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Methylation_Clustering/JHU_USC__HumanMethylation450.output_matrix.filtered_lambda0.2.CpGIslandSites.intersectingIDs.SomMat_IlluminaGA.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Douglas_LabMeetingImages/Meth450k_CIMPMarkerGenesOnly.intersectingIDs.tsv")