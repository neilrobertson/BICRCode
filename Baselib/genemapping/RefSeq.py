import csv
import sys
import os

class GenesMapping(dict):
    gfilesettings = {"id" : 0, "enstrans" : 1, "delimiter" : '\t' }
    
    def __init__(self, geneMappingFile):
        self.genes = self #dictionary of geneid -> [chr,strand,start,end]
        self.getGeneLocations(geneMappingFile)

    # read in the gene id mapping
    def getGeneLocations(self, geneidsFile):
        geneIdsInputFile = csv.reader(open(geneidsFile, "r"), delimiter=self.gfilesettings["delimiter"])
        for row in geneIdsInputFile:
            
            if row[0].startswith("#"):
                continue # skip header / comments
            
            geneid = row[self.gfilesettings["id"]]
            transcript = row[self.gfilesettings["enstrans"]]
            
            self[geneid] = transcript
        print "Genes locations filed:", len(self.genes)
