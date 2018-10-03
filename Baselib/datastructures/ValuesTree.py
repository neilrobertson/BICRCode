'''
Created on 13 Oct 2014

@author: johncole
'''

import time
import gc
import csv
from datastructures.genomeintervaltree import GenomeIntervalTree

def reverseOrder(line):
    (chrm,coord,value) = line
    return (chrm,coord,value)

def parseMethLine(line,reverse):
    if reverse:
        line = reverseOrder(line)
    
    (chrm,coord,value) = line
    chrm = chrm.strip()
    if not chrm.startswith("chr"):
        chrm = "chr"+chrm
    # coord should be in base 0 already, just need it as an int
    coord = int(coord)
    value = float(value)    
    return (chrm,coord,value)

class ValuesTree(GenomeIntervalTree):

    def __init__(self,methfilename,reverse=False):
        super(ValuesTree, self).__init__()        
        
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
            (chrm,coord,value) = parseMethLine(line,reverse)
            
            #Tests for a new chromosome
            if currentChrm != chrm:
                print chrm
                currentChrm = chrm
                #linesInserted = 0
            
            #Inserts a new node in the tree
            if not linesInserted == -1:
                #self.insertInterval(chrm, coord, coord+1, None)
                self.insertInterval(chrm, coord, coord+1, (value))
                linesInserted+=1
            
            #Prints to check whats going on
            if linesInserted % 100000 == 0:
                t2 = time.time()
                print "Inserted:",linesInserted,str((t2-t1)*1000.0) + ' ms'
                t1 = time.time()

        gc.enable()
            