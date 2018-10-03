import csv

class Pointlist(list):
    gfilesettings = {"chr" : 0, "point" : 1, "delimiter" : '\t' }
    
    def __init__(self, location):
        self.getPointlist(location)
    
    def getPointlist(self, location):
        inputFile = csv.reader(open(location, "r"), delimiter=self.gfilesettings["delimiter"])
        for row in inputFile:            
            if len(row)==0 or row[0].startswith("#"):
                continue # skip header / comments
            
            chr = row[self.gfilesettings["chr"]]
            
            if not chr.startswith("chr"):
                chr = "chr" + chr
            
            point = row[self.gfilesettings["point"]]
            
            self.append((chr, int(point)))
        print "Points locations filed:", len(self)

class PointlistWithValue(list):
    gfilesettings = {"chr" : 0, "point" : 1, "value" : 2, "delimiter" : '\t' }
    
    def __init__(self, location):
        self.getPointlist(location)
    
    def getPointlist(self, location):
        inputFile = csv.reader(open(location, "r"), delimiter=self.gfilesettings["delimiter"])
        for row in inputFile:            
            if row[0].startswith("#"):
                continue # skip header / comments
            
            chr = row[self.gfilesettings["chr"]]
            
            if not chr.startswith("chr"):
                chr = "chr" + chr
            
            point = row[self.gfilesettings["point"]]
            value = row[self.gfilesettings["value"]]
            
            self.append((chr, int(point),  float(value)))
        print "Points locations filed:", len(self)
