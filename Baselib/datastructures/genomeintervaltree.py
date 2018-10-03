from datastructures.cintervaltree import *
from collections import defaultdict
import sys

def intervalsOverlapOrAdjacent(astart,aend,bstart,bend):
    return aend >= bstart and astart <= bend

def intervalsOverlap(astart,aend,bstart,bend):
    return aend > bstart and astart < bend

class GenomeIntervalTree(object):
    def __init__(self):
        self.values = defaultdict(IntervalTree)
        self._entries = 0
    
    def status(self):
        for chrm in self.values:
            print chrm,":",self.values[chrm]
    
    #Does the recusion
    def insertInterval(self, chrm, start, stop, value):
        self.values[chrm].insert_interval(Interval(start, stop, value=value))
        self._entries += 1
    
    def hasChromosome(self, chrm):
        if not chrm.startswith("chr"):
            chrm = "chr" + chrm
        return chrm in self.values
    
    #Gets a tree for a given chromosome
    def getChromosomeIntervalTree(self, chrm):
        if not chrm.startswith("chr"):
            chrm = "chr" + chrm
        return self.values[chrm]

    #
    def getIntervalsInRange(self, chrm, start, stop):
        if not chrm.startswith("chr"):
            chrm = "chr" + chrm
        
        # get the interval tree for here
        if not self.hasChromosome(chrm):
            return [] # its on a chromosone we dont have any values for
        tree = self.getChromosomeIntervalTree(chrm)
        
        return tree.find_ucsc(min(start, stop), max(start, stop))
    
    
    def getValues(self,intervals):
        values = []
        for interval in intervals:
            values.append(interval.value)
            
        assert len(intervals)==len(values), "Values not same length as intervals"
        return values
    
    def getValuesInRange(self, chrm, start, stop):
        if not chrm.startswith("chr"):
            chrm = "chr" + chrm          
        
        return self.getValues(self.getIntervalsInRange(chrm, start, stop))

    def before(self, chrm, point, max_dist=None):
        # get the interval tree for here
        if not self.hasChromosome(chrm):
            return [] # its on a chromosone we dont have any values for
        tree = self.getChromosomeIntervalTree(chrm)

        if max_dist==None:
            max_dist = sys.maxint

        return tree.before(point+1, max_dist=max_dist)

    def after(self, chrm, point, max_dist=None):
        # get the interval tree for here
        if not self.hasChromosome(chrm):
            return [] # its on a chromosone we dont have any values for
        tree = self.getChromosomeIntervalTree(chrm)

        if max_dist==None:
            max_dist = sys.maxint

        return tree.after(point-1, max_dist=max_dist)        

    def __len__(self):
        return self._entries
