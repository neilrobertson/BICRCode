import csv
import os
import sys
import collections

netaffxBase = "~/mount/publicdata"

def __splitEntries__(string):
    return string.split("///")

def __splitParts__(string):
    return string.split("//")

def __splitString__(string):
    entries = []
    for entry in __splitEntries__(string):
        parts = []
        splitParts = __splitParts__(entry.strip())
        if len(splitParts)>1:
            for part in splitParts:
                parts.append(part.strip())
        else:
            parts = entry.strip()
        entries.append(parts)

    return entries

# Valid columns should be:
# ---
#Probe Set ID <--- This is the index and therefore isn't stored as a field
#GeneChip Array
#Species Scientific Name
#Annotation Date
#Sequence Type
#Sequence Source    
#Transcript ID(Array Design)
#Target Description
#Representative Public ID
#Archival UniGene Cluster
#UniGene ID
#Genome Version
#Alignments
#Gene Title
#Gene Symbol
#Chromosomal Location
#Unigene Cluster Type
#Ensembl
#Entrez Gene
#SwissProt
#EC
#OMIM
#RefSeq Protein ID
#RefSeq Transcript ID
#FlyBase   
#AGI
#WormBase
#MGI Name
#RGD Name
#SGD accession number
#Gene Ontology Biological Process
#Gene Ontology Cellular Component
#Gene Ontology Molecular Function
#Pathway
#InterPro
#Trans Membrane
#QTL
#Annotation Description
#Annotation Transcript Cluster
#Transcript Assignments
#Annotation Notes


cdfNames = {
            
            "Mouse430_2" : {
                            "default":"mm9",
                            "mm9": { "version" : "30", "type":"annot"}
                            }
            ,
            "HG-U133_Plus_2" : {
                               "default":"hg18",
                               "hg18": { "version" : "29", "type" : "annot"},
                               "hg19": {  "version" : "33", "type" : "annot"}
                               }
            }


class NetAffxAnnotation():
    def __init__(self, genome = None, cdfname=None):
        
        if cdfname == None:
            # assume HG-U133_Plus_2 array
            cdfname = "HG-U133_Plus_2"
        
        cdfname = cdfname.strip()
        assert cdfname in cdfNames, "'"+cdfname+"'"
        
        print cdfname
        print genome

        
        if genome == None:
            genome = cdfNames[cdfname]["default"]
        
        assert genome in cdfNames[cdfname]
        assert "version" in cdfNames[cdfname][genome]
        assert "type" in cdfNames[cdfname][genome]
    
        version = cdfNames[cdfname][genome]["version"]
        annottype = cdfNames[cdfname][genome]["type"]
        
        annotationFile = cdfname + ".na"+version+"."+annottype+".csv"
        fileLocation = os.path.expanduser(netaffxBase + "/" + genome + "/mappings/affy/" + annotationFile + "/" + annotationFile)
        
        inputFile = csv.reader(open(fileLocation, "r"), delimiter=',',  quotechar='"')
        
        self.symbolmap = None
        self.ensemblmap = None
        
        header = []
        self.meta = {}
        self.data = {} # mapping of header field -> {}
        self.keys = set()
        
        for row in inputFile:
            if row[0].startswith("#"):
                # comment or version info
                if row[0].startswith("#%"):
                    # this is meta annotation data
                    key, value = row[0][2:].split("=", 1)
                    self.meta[key] = value
                continue # done with this row - it was either a comment or we've stored the meta data
            
            if len(header) == 0:
                # we havent processed the header yet
                header = row
                
                for column in row:
                    self.data[column] = {} # this is the mapping that will store the data for each actual row
                
                continue # done with the header row
            
            # now we are onto actual data
            self.keys.add(row[0])
            for index in range(1, len(row)):
                value = row[index]
                
                self.data[header[index]][row[0]] = row[index]
            
        

    def getKeys(self):
        return self.keys

    def getMeta(self, meta):
        return self.meta[meta]

    def getValue(self, affy, column):
        return self.data[column][affy]

    def getValues(self, affy, column):
        
        if self.getValue(affy,column) == "---":
            return []
        else:
            return __splitString__(self.getValue(affy, column))
    
    def buildEnsemblMap(self):
        self.ensemblmap = collections.defaultdict(set)
        for affy in self.getKeys():
            for ensembl in self.getValues(affy, "Ensembl"):
                if ensembl != "---":
                    self.ensemblmap[ensembl].add(affy)
    
    def getAffysForEnsembl(self, ensemblId):
        if self.ensemblmap == None:
            self.buildEnsemblMap()
        
        if ensemblId in self.ensemblmap:
            return self.ensemblmap[ensemblId]
        else:
            return set()
    
    def getUniqueAffys(self,affys):
        uniqueChangingAffys = set()
        seenChangingGenes = set()
        for affy in affys:            
            genesymbol = self.getValue(affy, "Gene Symbol")            
            if genesymbol not in seenChangingGenes: 
                uniqueChangingAffys.add(affy)
                seenChangingGenes.add(genesymbol)        
        return uniqueChangingAffys 
    
    def buildSymbolMap(self):
        self.symbolmap = collections.defaultdict(set)
        for affy in self.getKeys():
            for affysymbol in self.getValues(affy, "Gene Symbol"):
                self.symbolmap[affysymbol].add(affy)
                self.symbolmap[affysymbol.upper()].add(affy)
    
    def getAffysForSymbol(self, symbol):
        if self.symbolmap == None:
            self.buildSymbolMap()
        
        if symbol in self.symbolmap:
            return self.symbolmap[symbol]
        elif symbol.upper() in self.symbolmap:
            return self.symbolmap[symbol.upper()]
        else:
            return set()

if __name__ == "__main__":
    annotation = NetAffxAnnotation()
    
    print annotation.getMeta("genome-version")
    
    print annotation.getValue("1007_s_at", "Gene Symbol")
    print annotation.getValue("1007_s_at", "Alignments")
    print annotation.getValue("1552411_at" ,"Gene Symbol")
    print annotation.getValue("207611_at" ,"Gene Symbol")
    
    print annotation.getValues("1007_s_at", "Gene Symbol")
    print annotation.getValues("1007_s_at", "Alignments")
    print annotation.getValues("1552411_at" ,"Gene Symbol")
