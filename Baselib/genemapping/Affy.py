import os
import sys

import csv
import genemapping.Ensembl as Ensembl

class GenesMapping(dict):
    
    mappings = {"hg18toU133P2" : "~/mount/publicdata/hg18/mappings/affyCoords/knownGenesOld.hg18toU133P2affy.csv"}
    
    gfilesettings = {"id" : 4, "chr" : 0, "start" : 2, "end":3, "strand":1, "delimiter" : '\t' }
    # chr1  -   4268    6628    225035_x_at
    
    def __init__(self, geneMappingFile = "hg18toU133P2"):
        self.genes = self # dictionary of geneid -> [chr,strand,start,end]
        assert geneMappingFile in self.mappings, geneMappingFile
        self.getGeneLocations(geneMappingFile)
    
    def addGene(self, geneid, chr, strand, start, end):
        self[geneid] = (chr, strand, start, end)
    
        # read in the gene id mapping
    def getGeneLocations(self, geneidsFile):
        geneIdsInputFile = csv.reader(open(os.path.expanduser(self.mappings[geneidsFile]), "r"), delimiter=self.gfilesettings["delimiter"])
        for row in geneIdsInputFile:
            
            if row[0].startswith("#"):
                continue # skip header / comments
            
            geneid = row[self.gfilesettings["id"]]
            chr = row[self.gfilesettings["chr"]]
            strand = row[self.gfilesettings["strand"]]
            start = row[self.gfilesettings["start"]]
            end = row[self.gfilesettings["end"]]
            
            self.addGene(geneid, chr, strand, int(start), int(end))

class FriendlyGeneNames(dict):
    def __init__(self, friendlyGenesFile):
        self.genes = self
        self.reverse = {}
        self.loadGenesToUse(friendlyGenesFile)

    def loadGenesToUse(self, friendlyGenesFile):
        genesToUseInputFile = csv.reader(open(friendlyGenesFile, "r"), delimiter='\t')
        for row in genesToUseInputFile:
            self[row[0]]=row[1]
            if row[1] not in self.reverse:
               self.reverse[row[1]]  = []
            self.reverse[row[1]].append(row[0])

    def getAffyIdsForGeneName(self, genename):
        if genename in self.reverse:
            return self.reverse[genename]
        else:
            return []
