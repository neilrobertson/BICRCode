'''
Created on 20 Jan 2012 

@author: mcbryan
'''
import time
import gc
import csv
from datastructures.genomeintervaltree import GenomeIntervalTree

def reverseOrder(line):
    (chrm,coord,lmeth,lunmeth,rmeth,runmeth) = line
    return (chrm,coord,rmeth,runmeth,lmeth,lunmeth)

def parseMethLine(line,reverse):
    if reverse:
        line = reverseOrder(line)
    
    (chrm,coord,lmeth,lunmeth,rmeth,runmeth) = line
    chrm = chrm.strip()
    if not chrm.startswith("chr"):
        chrm = "chr"+chrm
    # coord should be in base 0 already, just need it as an int
    coord = int(coord)
    lmeth=int(lmeth)
    lunmeth=int(lunmeth)
    rmeth=int(rmeth)
    runmeth=int(runmeth)        
    return (chrm,coord,lmeth,lunmeth,rmeth,runmeth)

class AggregateTree(GenomeIntervalTree):

    def __init__(self,methfilename,reverse=False):
        super(AggregateTree, self).__init__()        
        
        self.filename = methfilename
        self.loadValues(methfilename,reverse)
        
    def loadValues(self,filename,reverse):
        linesInserted = 0
        
        gc.disable()
        
        inputData = csv.reader(open(filename, "r"), delimiter='\t')
        
        currentChrm = None
        
        #Gets the time
        t1 = time.time()
        for line in inputData:
            if len(line)==0:
                continue
            
            #Parses each line
            (chrm,coord,lmeth,lunmeth,rmeth,runmeth) = parseMethLine(line,reverse)
            
            #Tests for a new chromosome
            if currentChrm != chrm:
                print chrm
                currentChrm = chrm
                #linesInserted = 0
            
            #Inserts a new node in the tree
            if not linesInserted == -1:
                #self.insertInterval(chrm, coord, coord+1, None)
                self.insertInterval(chrm, coord, coord+1, (lmeth,lunmeth,rmeth,runmeth))
                linesInserted+=1
            
            #Prints to check whats going on
            if linesInserted % 100000 == 0:
                t2 = time.time()
                print "Inserted:",linesInserted,str((t2-t1)*1000.0) + ' ms'
                t1 = time.time()

        gc.enable()
            