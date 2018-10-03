'''
Created on 2 Sep 2014

@author: neilrobertson
'''
from os import listdir
from os.path import isfile, join


class CNV_Low_Pass_DNASeq(object):
    def __init__(self, directory, outputDir, cancerGenomeAtlas, delimiter = "\t"):
        from CancerGenomeAtlas import CancerGenomeAtlas

        self.directory = CancerGenomeAtlas.checkDirectory(directory)
        assert self.directory
        self.cancerGenomeAtlas = cancerGenomeAtlas
        
        self.cnvFiles = [f for f in listdir(self.directory) if isfile(join(self.directory,f))]
        #print self.cnvFiles
        self.outputDelimiter = delimiter
        self.outputDir = outputDir
        
        
    def createSegmentationFile(self, outputFilename, arrayListFilename, genomeBuild = "hg19"):
        
        filteredFileList = [x for x in self.cnvFiles if x.find(".hg19.seg.txt") != -1]
        print len(filteredFileList)
        print "Creating segmentation files."
        includedSamples = []
        counter = 0
        with open(outputFilename, "w") as outputFile:
            
            for segFile in filteredFileList:
                sampleID = self.cancerGenomeAtlas.getSampleCodeFromFile(segFile).split("/")[0]
                print segFile
                print sampleID
                if sampleID not in includedSamples:
                    if sampleID.split("-")[-1] in ("01", "02", "03", "04", "05", "06", "07", "08", "09"):
                        segFilename = self.directory + r"/" + segFile
                        with open(segFilename, "r") as currentSegFile:
                            for i,row in enumerate(currentSegFile):
                                if i == 0: pass #header
                                else:
                                    outputFile.write(sampleID + self.outputDelimiter + row.strip() + "\n")
                        counter += 1
                        includedSamples.append(sampleID)
                        print "Completed writing segmentation file %s to global segmentation list. Filename: %s" % (str(counter), segFile)
        
        
        with open(arrayListFilename, "w") as sampleListFile:
            for includedSample in includedSamples:
                sampleListFile.write(includedSample + "\n")
                
        return outputFilename, arrayListFilename
    
    def dispose(self):    
        pass
        
    @staticmethod
    def getCNVDataFolderTypes():
        return ["HMS__IlluminaHiSeq_DNASeqC"]