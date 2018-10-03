'''
Created on 5 Sep 2014

@author: neilrobertson
'''

class ClusteredFile(object):
    def __init__(self):
        pass
    
    @staticmethod
    def mergeFiles(inputOne, inputTwo, outputFilename, annotationColumnKeys = 1, delimiter = "\t", fillMissing = None):
        from collections import defaultdict

        ids = defaultdict(list)
        dataLines = defaultdict(list)

        headerAnnotator = None
        
        firstContents = []
        inputDataLen_1 = None
        with open(inputOne, "r") as inputFileOne:
            for i, line in enumerate(inputFileOne):
                if i == 0: 
                    annotationKeys = delimiter.join(line.strip().split(delimiter)[:annotationColumnKeys])
                    sampleIDs = line.strip().split(delimiter)[annotationColumnKeys:]
                    headerAnnotator = annotationKeys
                    ids[annotationKeys].append(sampleIDs)
                    inputDataLen_1 = len(sampleIDs)
                else:
                    annotationKeys = delimiter.join(line.strip().split(delimiter)[:annotationColumnKeys])
                    dataPoints = line.strip().split(delimiter)[annotationColumnKeys:]
                    dataLines[annotationKeys].append(dataPoints)
                    firstContents.append(annotationKeys)
        
        secondContents = []
        inputDataLen_2 = None
        with open(inputTwo, "r") as inputFileTwo:
            for i, line in enumerate(inputFileTwo):
                if i == 0: 
                    annotationKeys = delimiter.join(line.strip().split(delimiter)[:annotationColumnKeys])
                    sampleIDs = line.strip().split(delimiter)[annotationColumnKeys:]
                    ids[annotationKeys].append(sampleIDs)
                    inputDataLen_2 = len(sampleIDs)
                else:
                    annotationKeys = delimiter.join(line.strip().split(delimiter)[:annotationColumnKeys])
                    dataPoints = line.strip().split(delimiter)[annotationColumnKeys:]
                    dataLines[annotationKeys].append(dataPoints)
                    secondContents.append(annotationKeys)
                        
        coexistingKeys = [x for x in firstContents if x in secondContents]
        nonCoexistingKeys_first = [x for x in firstContents if x not in secondContents]
        nonCoexistingKeys_second = [x for x in secondContents if x not in firstContents]
        print "We have %s co-existing keys!" % str(len(coexistingKeys))
        print "We have %s keys in the first file not in the second." % str(len(nonCoexistingKeys_first))
        print "We have %s keys in the second file not in the first." % str(len(nonCoexistingKeys_second))
        
        with open(outputFilename, "w") as output:
            
            header = headerAnnotator + delimiter
            for sampleIds in ids[headerAnnotator]:
                header += delimiter.join(sampleIds) + delimiter
            
            output.write(header.strip() + "\n")
            
            for dataPoint in coexistingKeys:
                dataLine = dataPoint + delimiter
                for dataPoints in  dataLines[dataPoint]:
                    dataLine += delimiter.join(dataPoints) + delimiter
                output.write(dataLine.strip() + "\n")
                
            if fillMissing:
                for dataPoint in nonCoexistingKeys_first:
                    dataLine = [dataPoint] + dataLines[dataPoint][0] + [0]*inputDataLen_2
                    output.write(delimiter.join(map(str, dataLine)) + "\n")
                for dataPoint in nonCoexistingKeys_second:
                    dataLine = [dataPoint] + [0]*inputDataLen_1 + dataLines[dataPoint][0]
                    output.write(delimiter.join(map(str, dataLine)) + "\n")   
        print "Completed file merge!"
        return outputFilename
               


    @staticmethod
    def buildClusterGroupDict(clusterAnnotationFile, delimiter = "\t"):
        assert clusterAnnotationFile
        
        annotationFile = open(clusterAnnotationFile, 'r')
        
        clusterDict = {}
        for i,row in enumerate(annotationFile):
            if i > 0:
                rowParts = row.strip().split(delimiter)
                annotationID = rowParts[0]
                clusterNo = rowParts[1]
                clusterDict[annotationID] = clusterNo
        return clusterDict
    
    @staticmethod
    def realignClusteredSampleIDs(inputFileName, outputFileName, clusteredOrderFilename, delimiter = "\t"):
        '''
        Realigns an input file so that the sample ids are in the order specified by the clusteredOrderFilename file - 
        a list of differently ordered Ids most likely the output of clustering
        '''
        from Utility.Utility import Utility
        
        assert inputFileName
        assert clusteredOrderFilename
        
        with open(inputFileName, 'r') as inputFile:
            
            idOrderFile = open(clusteredOrderFilename, 'r')
            idOrderRow = None
            for i, row in enumerate(idOrderFile):
                if i == 0:
                    idOrderRow = row
                else:
                    break
            orderedIds = idOrderRow.split(delimiter)
            orderedIds = tuple([x.strip() for x in orderedIds if not None]) 
            idOrderFile.flush()
            idOrderFile.close()
            
            if orderedIds:
                columns, indexToName = Utility.getColumns(inputFile)
                #nameToIndex = {v: k for k, v in indexToName.items()}
                failCount = 0
                nonIncludedSampledIds = []
                orderDict = {}
                for i in orderedIds:
                    if i not in indexToName.values():
                        failCount += 1
                        orderDict[i] = False
                        nonIncludedSampledIds.append(i)
                    else:
                        orderDict[i] = True
                print "There are %s samples in the input data set that exist in the output mutations data. %s non-included samples" % (str(len(orderedIds) - failCount) , str(failCount))
                
                outputFile = open(outputFileName, 'w') 

                headerLine = indexToName[0].strip() + delimiter + delimiter.join(orderedIds).strip() + "\n"
                outputFile.write(headerLine)
                lineCount = 0
                for geneSymbol in columns[indexToName[0]]:
                    line = geneSymbol + delimiter
                    for consecutiveId in orderedIds:
                        if orderDict[consecutiveId] == False:
                            line += str(-1) + delimiter
                        elif orderDict[consecutiveId] == True:
                            line += str(columns[consecutiveId][lineCount]) + delimiter
                    outputFile.write(line.rstrip(delimiter) + "\n")
                    lineCount += 1
                
                outputFile.flush()
                outputFile.close()
                print "Ordered file write completed..."
                print nonIncludedSampledIds
            else:
                print "The list of ordered IDs from the clustered file was not found... Shutting down.. beep..."
        return outputFileName 
    
#ClusteredFile.realignClusteredSampleIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/DNMT3B-CNV_Meth-Correlation/Douglas.CRC.copyNumber.DNMT3B.Meth450kSharedSamples.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/DNMT3B-CNV_Meth-Correlation/Douglas.CRC.copyNumber.DNMT3B.Meth450kSharedSamples.SamplesOrdered.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/PreProcessing/JHU_USC__HumanMethylation450.NA-Removed.SNPsRemoved.TumourSamples.MultiPlatformSharedSamples.Filtered-StdDev0.2.tsv", delimiter = "\t")
#ClusteredFile.realignClusteredSampleIDs("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/Clustering/MutationsHeatmaps/CRCCommonMutations.mutatedGeneBySample.MethClusterSharedSamples.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/Clustering/MutationsHeatmaps/CRCCommonMutations.mutatedGeneBySample.MethClusterSharedSamples.OrderedByCluster.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNA_Methylation/Clustering/MutationsHeatmaps/HumanMethylation450k.NoNoise.TumourSample.ClusteredSampleOrder.csv", delimiter = "\t")
#ClusteredFile.mergeFiles("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_11_02_15/Outputs/BCM__IlluminaGA_DNASeq.mutationSignatureBySample.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_11_02_15/Outputs/BCM__SOLiD_DNASeq.mutationSignatureBySample.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/SomaticMutations/AllSomaticMutations.mutationSignatureBySample.tsv", annotationColumnKeys = 10, delimiter = "\t", fillMissing = True)
#ClusteredFile.mergeFiles("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_11_02_15/Outputs/BCM__IlluminaGA_DNASeq.mutatedGeneBySample.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_11_02_15/Outputs/BCM__SOLiD_DNASeq.mutatedGeneBySample.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/SomaticMutations/AllSomaticMutations.mutatedGeneBySample.tsv", annotationColumnKeys = 1, delimiter = "\t", fillMissing = True)
       
#ClusteredFile.realignClusteredSampleIDs("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Methylation_BoxPlots/JHU_USC__HumanMethylation450.output_matrix.IDsOnly.intersectingIDs.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Methylation_BoxPlots/JHU_USC__HumanMethylation450.output_matrix.IDsOnly.intersectingIDs.alignedToCluster.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Methylation_Clustering/JHU_USC__HumanMethylation450.output_matrix.filtered_lambda0.2.CpGIslandSites.intersectingIDs.SomMat_IlluminaGA.tsv.hclust.cluster.tsv")            
#ClusteredFile.realignClusteredSampleIDs("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Douglas_LabMeetingImages/Meth450k_CIMPMarkerGenesOnly.intersectingIDs.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Douglas_LabMeetingImages/Meth450k_CIMPMarkerGenesOnly.intersectingIDs.alignedToCluster.tsv", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Methylation_Clustering/JHU_USC__HumanMethylation450.output_matrix.filtered_lambda0.2.CpGIslandSites.intersectingIDs.SomMat_IlluminaGA.tsv.hclust.cluster.tsv")        