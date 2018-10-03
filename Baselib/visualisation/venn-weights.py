import math
import numpy
import csv
import sys
import getopt
import gc
import string
import collections

import numpy as np


gc.disable()

from bed.treatment import SimpleBed
from bed.treatment import Bed as BedIntervalTree
from datastructures.genomeintervaltree import GenomeIntervalTree
from itertools import chain, combinations

def getVenneulerString(vennWeights):
    assert isinstance(vennWeights, object)
    return """
library(eulerr)

v <- euler(c(%s))
v$labels <- c(<INSERT BED LABELS>)
    
png("<INSERT FILENAME>",width=960,height=960,pointsize=24)
plot(v)
dev.off()
""" % vennWeights

class Interval(object):
    def __init__(self,chrm,start,end):
        self.chrm = chrm
        self.start = start
        self.end = end

def mergeIntervals(intervals):
        sortedList = [Interval(chrm,start,end) for chrm,start,end in intervals] 
        
        #Sorts the list
        sortedList.sort(key=lambda x: (x.chrm, x.start, x.end))
        
        mergedIntervals = []
         
        #merges the list
        previous = None
        for current in sortedList:
            
            if previous == None:
                previous = current
                continue

            elif (current.end >= previous.start and current.start <= previous.end and previous.chrm == current.chrm):
                
                start = min(current.start, previous.start)
                end = max(current.end, previous.end)
                
                previous = Interval(previous.chrm,start,end)
               
            else:
                mergedIntervals.append(previous)
                previous = current
        
        mergedIntervals.append(previous)
        
        return mergedIntervals
    
def intervalListToFile(path,intervalList):
    
    f = open(path + "/temp.csv", "w")
    for i in intervalList:
        f.write(str(str(i.chrm) + "\t" + str(i.start) + "\t" + str(i.end) + "\n"))
    f.close
    
    return path + "/temp.csv"  
    
def determineOverlap(regions,treatment):
   
    overlappingIntervals = []    
    overlappingBP = 0
    
    for row in regions:
        (chr, start, stop) = row
        intervals = treatment.getIntervalsInRange(chr, start, stop)
            
        for interval in intervals:
            # how much of interval is overlapped with chr, start, stop
            overlappingBP += min(interval.end,stop) - max(interval.start,start)
            
            overlappingIntervals.append((chr,max(interval.start,start),min(interval.end,stop)))
            
    return overlappingBP, overlappingIntervals

def numbBP(regions):
    bp = 0
    for chr,start,stop in regions:
        bp += stop-start
    return bp


def compareTwoFiles(a,b):
    
    aBed = SimpleBed(a)
    bTree = BedIntervalTree(b)

    return determineOverlap(aBed,bTree)
    
def compareListAndFile(ab,c):
    
    cTree = BedIntervalTree(c)
    
    return determineOverlap(ab,cTree)


def concatenateSimpleBed (fileList):
    
    a = []
    for i in fileList:
        a.extend(SimpleBed(i))
        
    return mergeIntervals(a)

class VennDiagram:

    def __init__(self, name_set_tuples_or_dict):
        if hasattr(name_set_tuples_or_dict, 'items'):
            sets = name_set_tuples_or_dict.items()
        else:
            sets = name_set_tuples_or_dict
        self.sets = []
        for name, group in sets:
            self.sets.append((name, group))

    default_colors = [ (228 / 255.0, 26 / 255.0, 28 / 255.0),
                              (55 / 255.0, 126 / 255.0, 184 / 255.0),
                              (77 / 255.0, 175 / 255.0, 74 / 255.0)]

    def plot_normal(self, output_filename, width=8):
        return self._venn_plot_weights(output_filename, width, width)

    def _get_set_names(self):
        return [name for name,values in self.sets]

    def _get_set_values(self):
        return [values for name,values in self.sets]
    

    def _venn_plot_weights(self, output_filename, width=8, height=8):
        """Plot a venn diagram into the pdf file output_filename.
        Takes a dictionary of sets and does the intersection calculation in python
        (which hopefully is a lot faster than passing 10k set elements to R)
        """
        weights = []
        intersects = []
        excludes = []
        sets_by_power_of_two = {}
            
        for ii, (kset, iset) in enumerate(self.sets):
            sets_by_power_of_two[2**ii] = iset
            
        for i in xrange(1, 2**len(self.sets)):
            sets_to_intersect = []
            to_exclude = []
            for ii in xrange(0, len(self.sets)):
                if (i & (2**ii)):
                    sets_to_intersect.append(sets_by_power_of_two[i & (2**ii)])
                else:
                    to_exclude.append(sets_by_power_of_two[(2**ii)])
            intersects.append(sets_to_intersect)
            excludes.append(to_exclude)
        
        return intersects, excludes

if __name__ == '__main__':
    
    import sys

    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:")
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err)
        sys.exit(2)

    # add executing directory as part of path
    sys.path.append(sys.path[0])
    
    sets = list()
    raweights = []
    weights = []
    counter = 0
    files = []
    
    for opt, arg in opts:
        if opt == "-s":
            sets.append((str(counter),arg))
            files.append(arg)
            counter +=1
    
    letters = list(string.ascii_uppercase)
    vennLabels = dict(zip(files, letters))
         
    #Gets the intersect and excluded lists
    vd = VennDiagram(sets)
    
    intersect,exclude = vd.plot_normal("", 1)
    
    for i in range(0,len(intersect)):
        
        intList = intersect[i]
        excList = exclude[i]
        
        print intList
        
        #Test for excluded regions:
        excludeUnionBed = None
        excludeBp = 0
        
        if (len(excList)) >0:
            #Concatenates the excluded bed list and gets the union of its regions
            excludeUnionBed = concatenateSimpleBed(excList)
            #Outputs the union list to a file
            excludeUnionBed = intervalListToFile(sys.path[0],excludeUnionBed)
        
        if len(intList) == 1:
            overlap = numbBP(SimpleBed(intList[0]))
            raweights.append(overlap)
            if excludeUnionBed != None:
                excludeBp,_ = compareTwoFiles(intList[0],excludeUnionBed)
            weights.append(overlap - excludeBp)
            
        if len(intList) == 2:
            overlap,bed = compareTwoFiles(intList[0],intList[1])
            raweights.append(overlap)
            if excludeUnionBed != None:
                excludeBp,_ = compareListAndFile(bed,excludeUnionBed)
            weights.append(overlap - excludeBp)

        if len(intList) > 2:
            _,bed = compareTwoFiles(intList[0],intList[1])
            
            for j in range(2,len(intList)):
                overlap,bed = compareListAndFile(bed,intList[j])
                
                if j == len(intList)-1:
                    raweights.append(overlap)
                    if excludeUnionBed != None:
                        excludeBp,_ = compareListAndFile(bed,excludeUnionBed)
                    weights.append(overlap - excludeBp)
        

    print "Raw Weights:"
    for i in raweights:
         print i  
    
    print "Weights:"
    for i in weights:
        print i
        
    vennerable = "0"   
    for i in weights:
        vennerable += ", " + str(i) 
    print vennerable
    
    # Writes overlap details to console in venneulur R code
    vennDetails = {}
    vennOrder = []
    for i in range(0,len(intersect)):
        intList = intersect[i]        
        intList = [vennLabels[x] for x in intList]
        labels = "&".join(intList)
        vennDetails[labels] = weights[i]
        vennOrder.append(labels)
    
    print files
    print "\n"
    venneulur = []
    for item in vennOrder:
        venneulur.append("\"%s\"=%s" % (item, vennDetails[item]))
    print getVenneulerString(",".join(venneulur))        
            
    
    
    
   
 