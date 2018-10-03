import csv
import sys

class GeneList(list):
    def __init__(self, genesToUseFileLocation):
        self.seengenes = set()
        self.genesToUseFileLocation = genesToUseFileLocation
        self.genes = self
        self.loadGenesToUse(genesToUseFileLocation)

    def loadGenesToUse(self, genesToUseFile):
        genesToUseInputFile = csv.reader(open(genesToUseFile, "r"), delimiter='\t')
        for row in genesToUseInputFile:
            # skip commented rows
            if len(row) == 0 or row[0].startswith("#") or row[0].startswith(">"):
                continue
            for column in row:
                if column.startswith("#") or column.startswith(">") or column == "":
                    continue
                for gene in column.split(","):
                    if not gene in self.seengenes and gene != "":
                        self.append(gene)
                        self.seengenes.add(gene)
        #print >> sys.stderr, "Genes To Use:", len(self)

    def getFriendlyName(self):
        return self.genesToUseFileLocation[self.genesToUseFileLocation.rfind("/")+1:]
    
    def getFullName(self):
        return self.genesToUseFileLocation

class GeneListWithValue(list):
    def __init__(self, genesToUseFileLocation):
        self.genesToUseFileLocation = genesToUseFileLocation
        self.loadGenesToUse(genesToUseFileLocation)

    def loadGenesToUse(self, genesToUseFile):
        genesToUseInputFile = csv.reader(open(genesToUseFile, "r"), delimiter='\t')
        for row in genesToUseInputFile:
            self.append((row[0], float(row[1])))
        #print "Genes To Use:", len(self)

    def getFriendlyName(self):
        return self.genesToUseFileLocation[self.genesToUseFileLocation.rfind("/")+1:]

    def getFullName(self):
        return self.genesToUseFileLocation
