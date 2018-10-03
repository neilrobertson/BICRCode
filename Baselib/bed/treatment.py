import csv
from datastructures.genomeintervaltree import GenomeIntervalTree
import gc
import os
import time
import sys

from visualisation.EasyDialogs import ProgressBar

# default strand is +
class BedEntry():
    def __init__(self,chr,start,stop,data,strand="+"):
        assert strand in ["+","-"]
        if not chr.startswith("chr"):
            chr = "chr" + chr
        self.chr = chr
        self.start = int(start)
        self.stop = int(stop)
        self.data = data
        self.strand = strand
        assert int(start)<=int(stop),  "Start should be before stop:"+str(start)+"-"+str(stop)
    
    def __str__(self):
        return self.chr +":"+str(self.start)+"-"+str(self.stop)

class BedFile():
    def __init__(self, fileName, chrPos = 0,  startPos = 1,  stopPos = 2, padding = 0, strandPos = 4, defaultkeys = []):
        self.fileName = fileName
        self.keys = defaultkeys
        self.chrPos = chrPos
        self.startPos = startPos
        self.stopPos = stopPos
        self.strandPos = strandPos
        self.padding = padding
        
    def __iter__(self):
        for row in csv.reader(open(self.fileName, "r"), delimiter='\t'):
            if len(row)==0:
                continue
            elif row[0].startswith("#"):
                for key in row:
                    self.keys.append(key.replace("#",""))
                print "Loading BED file with keys: "+repr(self.keys)
                continue
            #do int() casting when reading -> get it done early on

            chr = row[self.chrPos]
            
            if not chr.startswith("chr"):
                chr = "chr" + chr
            
            start = int(row[self.startPos])
            stop = int(row[self.stopPos])
            
            strand = "+"
            if len(row) >= self.strandPos + 1 and row[self.strandPos] in ["+","-"]:
                strand = row[self.strandPos]
            
            if len(self.keys) == 0:
                data = row
            else:
                data = {}
                for i in range(0, len(row)):
                    if row[i]=="":
                        continue
                    data[self.keys[i]]=row[i]
            
            assert start<=stop,  "Start should be before stop:"+str(start)+"-"+str(stop)
            yield BedEntry(chr, start, stop, data, strand = strand)
    
    def getFriendlyName(self):
        return self.fileName[self.fileName.rfind("/")+1:]
    
    def getFullName(self):
        return self.fileName
    
    def buildIntervalTree(self):
        intervalTree = GenomeIntervalTree()
        for entry in self:
            intervalTree.insertInterval(entry.chr, entry.start, entry.stop, entry.data)
        return intervalTree
        

class SimpleBed(list):
    def __init__(self, treatmentFileName,  padding = 0, validChrs = None):
        self.loadValues(treatmentFileName,  padding, validChrs)
    
    def loadValues(self,  valuesFile,  padding = 0, validChrs = None):
        count = 0
        valuesInputFile = csv.reader(open(valuesFile, "r"), delimiter='\t')
        self.header=["chr","start","stop"]
        self.extradata = {}
        for row in valuesInputFile:
            # skip empty rows
            if len(row)==0:
                continue
            # if row[0] starts with a # its a header row and defines our column names
            if row[0].startswith("#"):
                self.header = row
                continue
            #print row[0], row[1], row[2]
            chr = row[0]
            if not chr.startswith("chr"):
                chr = "chr" + chr
            
            if validChrs != None:
                if chr not in validChrs:
                    continue
            
            count = count + 1
            
            #do int() casting when reading -> get it done early on
            start = int(row[1]) - padding
            stop = int(row[2]) + padding
            assert start<=stop,  "Start should be before stop:"+str(start)+"-"+str(stop)
            
            data = (chr, start, stop)
            
            self.append(data)
            
            if len(row)>3:
                # we have extra data
                # we should have the same number of columns as we have in the header 
                if len(row)==len(self.header):
                    extradata = {}
                    for i in range(3,len(self.header)):
                        extradata[self.header[i]]=row[i]
                    assert data not in self.extradata, "Duplicates in SimpleBed not allowed " + str(row) 
                else:
                    extradata = row[3:]                   
                self.extradata[data]=extradata 
            
        #print "Bed Values filed:", count
        return self

    def getExtraData(self,data):
        return self.extradata[data] if data in self.extradata else {}
    
    def getTree(self):
        return SimpleBedTree(self)
    
class SimpleBedTree(GenomeIntervalTree):
    def __init__(self,simplebed,extend=None,padding = 0):
        for chrm,start,stop in simplebed:
            
            assert extend==None, "Not implemented yet"

            # insert padding as required
            start -= padding
            stop += padding
            
            assert start<=stop,  "Start should be before stop:"+str(start)+"-"+str(stop)
            
            self.insertInterval(chrm, start, stop, 1)
        
        

from subprocess import check_output

def wc(filename):
    return int(check_output(["wc", "-l", filename]).split()[0])

class Bed(GenomeIntervalTree):
    def __init__(self, treatmentFileName,limit=None,padding = 0,extend = None):        
        super(Bed, self).__init__()
        
        self.count = 0
        self.treatmentFileName = treatmentFileName
        self.loadValues(treatmentFileName,limit,padding,extend)
        
    def getFriendlyName(self):
        return self.treatmentFileName[self.treatmentFileName.rfind("/")+1:]
    
    def getFullName(self):
        return self.treatmentFileName

    # get the values from the treatment
    def loadValues(self, valuesFile, limit = None, padding = 0,extend = None,useProgressBar = ""):
 
        total = wc(valuesFile)
        
        if useProgressBar == None and total > 1000000:
            useProgressBar = True
            pb = ProgressBar("Loading "+self.getFriendlyName(),total)
        
        with open(valuesFile,"r") as f:
            # temporarily shut down the gc as we will be potentially doing lots of allocations
            # reenable once done loading
            gc.disable()
            valuesInputFile = csv.reader(f, delimiter='\t')
            
            try:
                for row in valuesInputFile:

                    if len(row)==0 or row[0].startswith("#"):
                        continue
                    
                    self.count = self.count + 1
                    
                    if useProgressBar:
                        if self.count % 100000 == 0:
                            pb.label(str(self.count) + " of " + str(total) + " processed ("+"{0:0.1f}".format(100.0*float(self.count)/float(total))+" %)")
                            pb.set(self.count)
                    
                    chrm = row[0]
                    if not chrm.startswith("chr"):
                        chrm = "chr" + chrm
                    
                    if len(row) < 2:
                        print row 
                        print valuesFile
                        print len(row)
                        
                    start = int(row[1])
                    stop = int(row[2])
                    
                    # apply strand specific extension to the read / feature
                    if extend != None:
                        assert len(row)>=6
                        strand = row[5]
                        assert strand in ['+','-']
                        
                        if strand == "+":
                            stop = start + extend
                        else:
                            start = stop - extend
                    
                    # insert padding as required
                    start -= padding
                    stop += padding
                    
                    assert start<=stop,  "Start should be before stop:"+chrm+str(start)+"-"+str(stop)
                    
                    self.insertInterval(chrm, start, stop, 1)
                    if limit != None and self.count > limit: # debugging quickness
                        break
            except csv.Error, e:
                print row
                sys.exit('file %s, line %d: %s' % (valuesFile, valuesInputFile.line_num, e))
                    
        if useProgressBar:
            pb.__del__()
        gc.enable()
        #print "Bed Values filed:", count

class Vstep(GenomeIntervalTree):
    def __init__(self, treatmentFileName,width=200,limit=0):        
        super(Vstep, self).__init__()
        self.count = 0
        self.width = width
        self.treatmentFileName = treatmentFileName
        self.loadValues(treatmentFileName,limit)
        
    # get the values from the treatment
    def loadValues(self, valuesFile, limit = 0):
        values = {} # dictionary of chr -> interval tree (start, end, value)
        valuesInputFile = csv.reader(open(valuesFile, "r"), delimiter='\t')
        for row in valuesInputFile:
            if len(row)==0 or row[0].startswith("#"):
                continue
            #do int() casting when reading -> get it done early on
            self.count = self.count + 1
            chr = row[0]
            
            if not chr.startswith("chr"):
                chr = "chr" + chr
            
            start = int(row[1])
            stop = start + self.width
            val = int(row[2])
            self.insertInterval(chr, start, stop, val)
            if limit > 0 and self.count > limit: # debugging quickness
                break
        #print "Bed Values filed:", count
        return values


# Assumes a line structure like this:
# chr1    1000    2000    value
# put the "value"s into the tree
# I'm not sure that's technically a bed, but it is useful.
# tries to force the value into a float, at the moment.
class ValueBed(GenomeIntervalTree):
    def __init__(self, treatmentFileName, padding=0):
        super(ValueBed, self).__init__()
        self.count = 0
        self.treatmentFileName = treatmentFileName
        self.loadValues(treatmentFileName, padding)
    
    # get the values from the treatment
    def loadValues(self, valuesFile, padding):
        values = {} # dictionary of chr -> interval tree (start, end, value)
        valuesInputFile = csv.reader(open(valuesFile, "r"), delimiter='\t')
        for row in valuesInputFile:
            if len(row)==0 or row[0].startswith("#"):
                continue
            self.count = self.count + 1
            chr = row[0]
            
            if not chr.startswith("chr"):
                chr = "chr" + chr
            
            start = int(row[1]) - padding
            stop = int(row[2]) + padding
            try:
                value = float(row[3])
            except ValueError:
                value = row[3]
            assert start<=stop,  "Start should be before stop:"+str(start)+"-"+str(stop)
            self.insertInterval(chr, start, stop, value)
        #print "Bed Values filed:", count
        return values

class ExtendedBed(GenomeIntervalTree):
    def __init__(self, treatmentFileName,chrPos = 0, startPos = 1, stopPos = 2,limit=0, defaultkeys = ["chr","start","stop"],forcekeys=False):
        super(ExtendedBed, self).__init__()
        self.loadValues(treatmentFileName,chrPos, startPos, stopPos, limit, defaultkeys, forcekeys)

    # get the values from the treatment
    def loadValues(self, valuesFile, chrPos = 0,  startPos = 1,  stopPos = 2, limit = 0, defaultkeys = ["chr","start","stop"],forcekeys=False):
        values = {} # dictionary of chr -> interval tree (start, end, value)
        count = 0
        self.keys = defaultkeys
        valuesInputFile = csv.reader(open(valuesFile, "r"), delimiter='\t')
        for row in valuesInputFile:
            if len(row)==0:
                continue
            elif row[0].startswith("#"):
                if not forcekeys:
                    self.keys = []
                    for key in row:
                        self.keys.append(key.replace("#",""))
                    print "Loading BED file with keys: "+repr(self.keys)
                continue
            #do int() casting when reading -> get it done early on
            count = count + 1
            chr = row[chrPos]
            
            if not chr.startswith("chr"):
                chr = "chr" + chr
            
            start = int(row[startPos])
            stop = int(row[stopPos])
            
            data = {}
            
            for i in range(0, len(row)):
                if row[i]=="":
                    continue
                data[self.keys[i]]=row[i]
            
            assert start<=stop,  "Start should be before stop:"+str(start)+"-"+str(stop)
            self.insertInterval(chr, start, stop, data)
            if limit > 0 and count > limit: # debugging quickness
                break
        
        print valuesFile[valuesFile.rfind("/")+1:]+ ": Bed Values filed:", count
        return values
