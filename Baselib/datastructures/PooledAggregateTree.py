'''
Created on 25 Jul 2014

@author: neilrobertson
'''
import time
import gc
import csv
from datastructures.genomeintervaltree import GenomeIntervalTree

def reverseOrder(line):
    #This needs correctly implemented
    (chrm,coord,lmeth,lunmeth,rmeth,runmeth) = line### INCORRECT
    return (chrm,coord,rmeth,runmeth,lmeth,lunmeth)### INCORRECT

def parseMethLine(line,reverse):
    if reverse:
        line = reverseOrder(line) # Not correctly implemented
    
    chrm = line[0].strip()
    if not chrm.startswith("chr"):
        chrm = "chr"+chrm
    # coord should be in base 0 already, just need it as an int
    coord = int(line[1])
    methylationValues = line[2:]
    return (chrm,coord,methylationValues)

class PooledAggregateTree(GenomeIntervalTree):

    def __init__(self,methfilename,reverse=False):
        super(PooledAggregateTree, self).__init__()        
        
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
            (chrm,coord,methylationValues) = parseMethLine(line,reverse)
            #Tests for a new chromosome
            if currentChrm != chrm:
                print chrm
                currentChrm = chrm
                #linesInserted = 0
            
            #Inserts a new node in the tree
            if not linesInserted == -1:
                #self.insertInterval(chrm, coord, coord+1, None)
                self.insertInterval(chrm, coord, coord+1, methylationValues)
                linesInserted+=1
            #Prints to check whats going on
            if linesInserted % 100000 == 0:
                t2 = time.time()
                print "Inserted:",linesInserted,str((t2-t1)*1000.0) + ' ms'
                t1 = time.time()

        gc.enable()
            