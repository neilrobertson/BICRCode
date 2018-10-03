'''
Created on 16 Sep 2014

@author: neilrobertson
'''
import os

class DataFiles(object):
    def __init__(self):
        pass

    @staticmethod
    def getSampleIDs_forFiltering():
        filename = "SampleIDsFiltering.txt"
        data = DataFiles.loadDataFile(filename)
        return data
    
    @staticmethod
    def getClinicalPatientIDs():
        filename = "ClinicalPatientIDs.txt"
        data = DataFiles.loadDataFile(filename)
        return data
    
    @staticmethod
    def getGenesList():
        filename = "GenesOfInterest.txt"
        data = DataFiles.loadDataFile(filename)
        data = [x.upper() for x in data]
        return data
    
    @staticmethod
    def getCpGReferences():
        filename = "CpGReferencesMatchToPosition.txt"
        data = DataFiles.loadDataFile(filename)
        return data
    
    @staticmethod
    def getCpGSitesOfInterest_Meth450k():
        filename = "MethylationArray_SitesOfInterest.txt"
        ACTION_LABEL = r"CpGSitesOfInterest"
        MAINTAIN_SET = True
        data = DataFiles.loadDataFile(filename)
        return data, MAINTAIN_SET, ACTION_LABEL
    
    @staticmethod
    def getCpGIslandSites_Meth450k():
        filename = "Methylation450kArray_CpGIslandSites.txt"
        ACTION_LABEL = r"CpGIslandSites"
        MAINTAIN_SET = True
        data = DataFiles.loadDataFile(filename)
        return data, MAINTAIN_SET, ACTION_LABEL
    
    
    @staticmethod
    def getSNPSites_Meth450k():
        filename = "Methylation450kArray_SNP_CpGs.txt"
        ACTION_LABEL = r"SNPsRemoved"
        MAINTAIN_SET = False
        data = DataFiles.loadDataFile(filename)
        return data, MAINTAIN_SET, ACTION_LABEL
    
    
    @staticmethod
    def loadDataFile(filename):
        data = []
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        with open(os.path.join(__location__, filename), 'r') as inputFile:
            for row in inputFile:
                if row.strip() not in ("", None):
                    data.append(row.strip())
        return data
    
    @staticmethod
    def openFile(filename):
        fileContent = ""
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        with open(os.path.join(__location__, filename), 'r') as inputFile:
            fileContent = inputFile.readlines()
        return fileContent