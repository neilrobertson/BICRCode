'''
Created on 25 Jul 2012

@author: mcbryan
'''

from vstep.treatment import Vstep as VstepFile
from bed.treatment import Bed as BedFile
from genemapping.chrmEnds import ChromosomeEnds

#####

# How Missing Values are dealt with - mostly relevant to vstep files
def missingValuesDontCount(chrm,values,start,end):
    if len(values)==0:
        return []
    else:
        return values
    
def missingValuesCountAsZero(chrm,values,start,end):
    if len(values)==0:
        return [0]
    else:
        return values

#####

class DataBehaviour(object):
    
    def __init__(self):
        pass
    
    # used to merge etc values together (raw value specific)
    def valuesBehaviour(self,chrm,values,start,end):
        assert False, "Values Behaviour not defined"

    # used to get values for a slice (raw value specific)
    def getValues(self,chrm,start,end):
        assert False, "Get Values not defined"
    
    # used for working out the bounds on heatmaps
    def heatmapHasNegativeValues(self):
        return False
    
    # used for working out the bounds on heatmaps
    def heatmapUpperLowerBounds(self):
        return None,None


#####
# VStep Files are pzs approach
#####
    
class Vstep(DataBehaviour):
    # vstep uses avg([...], [...], ...)
    # vstep-processed file (chr start end value)
    
    def __init__(self,filename,vstepWidth,zeroes = False):
        self.treatment = VstepFile(filename, 200 if vstepWidth == None else vstepWidth)
        
        self.getValues = self.treatment.getValuesInRange
        
        if zeroes:        
            self.valuesBehaviour = missingValuesCountAsZero
        else:
            self.valuesBehaviour = missingValuesDontCount

#####

class VstepIntermediateAverage(DataBehaviour):
    # vstep intermediate uses avg(avg([...]), avg([...]), ...)
    # vstep-processed file (chr start end value)
    
    def __init__(self,filename,vstepWidth,intermediateAvgFunc,zeroes=False):
        self.treatment = VstepFile(filename, 200 if vstepWidth == None else vstepWidth)
        self.intermediateAvgFunc = intermediateAvgFunc
        self.getValues = self.treatment.getValuesInRange

        if zeroes:
            self.missingValues = missingValuesCountAsZero
        else:
            self.missingValues = missingValuesDontCount

    def valuesBehaviour(self,chrm, values, start, end):
        if len(values)==0:
            return self.missingValues(chrm,values,start,end)
        else:
            return [self.intermediateAvgFunc(values)]
        
#####
        
import math

#For a BED file input
class Bed(DataBehaviour):
    
    #Constructor
    def __init__(self,filename,extend=None,fractions=True,normalizeLibrarySize=True):
        self.data = BedFile(filename,extend=extend)
        self.getValues = self.data.getIntervalsInRange
        self.fractions = fractions
        self.normalizeLibrarySize = normalizeLibrarySize # divide by 1mil to get a normalization factor
    
    def getLibraryNormalizationFactor(self):
        return self.data.count / 1000000
    
    # bed uses avg(  sum([...]) / bp , sum([...]) / bp, ...   )
    def fractionallyCountByBP(self,chrm, intervals, start, end):
        # sum = 0 if no values
        #assert math.fsum(intervals)>=0.0,  "Sum negative"
        assert end > start, str(end) + " > " + str(start)
        
        values = []
        for interval in intervals:
            assert interval.end>interval.start
            
            containedlength = min(interval.end,end) - max(interval.start,start)
            length = interval.end - interval.start
            
            values.append(float(containedlength) / float(length))
            
        # normalize by size of range
        value = math.fsum(values) / float(end-start)
        
        # normalize by library size
        if self.normalizeLibrarySize:
            value /= self.getLibraryNormalizationFactor()
        
        return value
    
    # avg( sum([...]), sum([...]))
    def count(self,chrm, values, start, end):
        value = len(values)        
        if self.normalizeLibrarySize:
            value /= self.getLibraryNormalizationFactor()        
        return value
    
    def valuesBehaviour(self,chrm,values,start,end):
        if self.fractions:
            return [self.fractionallyCountByBP(chrm,values,start,end)]
        else:
            return [self.count(chrm,values,start,end)]

#####

# extends Bed from above so we can use fractionallyCountByBP + count
class BedWithControl(Bed):
    
    # sample is a Bed from above
    def __init__(self,controlFilename,sample,extend=None,fractions=True,normalizeLibrarySize=True):
        
        # control is already loaded in the sample Bed object
        self.sample = sample
        self.control = BedFile(controlFilename,extend=extend)
        
        self.fractions = fractions
        self.normalizeLibrarySize = normalizeLibrarySize # divide by 1mil to get a normalization factor
        
        assert self.fractions == self.sample.fractions, "Bed and BedWithControl in different formats"
        assert self.normalizeLibrarySize == self.sample.normalizeLibrarySize
    
    def getValues(self,chrm, start, stop):
        # store a tuple of sample/cntrl values
        return self.sample.getValues(chrm,start,stop),self.control.getIntervalsInRange(chrm, start, stop)
    
    def getLibraryNormalizationFactor(self):
        return self.control.count / 1000000
    
    def valuesBehaviour(self,chrm,values,start,end):
        
        sampleValues, controlValues = values
        
        if self.fractions:
            return [self.sample.fractionallyCountByBP(chrm,sampleValues,start,end) - self.fractionallyCountByBP(chrm,controlValues,start,end)]
        else:
            return [self.sample.count(chrm,sampleValues,start,end) - self.count(chrm,controlValues,start,end)]
        
#####

#For a BED file input (where we only require a count)
class BedCount(DataBehaviour):
    
    #Constructor
    def __init__(self,filename,extend=None,fractions=True,normalizeLibrarySize=False):
        self.data = BedFile(filename,extend=extend)
        self.getValues = self.data.getIntervalsInRange
        self.fractions = fractions
    
    # bed uses avg(  sum([...]) / bp , sum([...]) / bp, ...   )
    def fractionallyCountByBP(self,chrm, intervals, start, end):
        # sum = 0 if no values
        #assert math.fsum(intervals)>=0.0,  "Sum negative"
        assert end > start, str(end) + " > " + str(start)
        
        values = []
        for interval in intervals:
            
            assert interval.end > interval.start
            
            containedlength = min(interval.end,end) - max(interval.start,start)
            
            values.append(float(containedlength))
            
        # normalize by size of range
        value = math.fsum(values) / float(end-start)
        
        #print value
        return value
    
    # avg( sum([...]), sum([...]))
    def count(self,chrm, values, start, end):
        value = len(values)              
        return value
    
    def valuesBehaviour(self,chrm,values,start,end):
        if self.fractions:
            return [self.fractionallyCountByBP(chrm,values,start,end)]
        else:
            return [self.count(chrm,values,start,end)]

#####


#For a BED file input (where we require the average of a series of values)
from datastructures.ValuesTree import ValuesTree

class BedValues(DataBehaviour):
    
    #Constructor
    def __init__(self,filename,extend=None,fractions=True,normalizeLibrarySize=False):
        self.treatment = ValuesTree(filename)
        self.valuesBehaviour = missingValuesDontCount

    #Gets the mean value
    def meanValue(self,values):
        
        total = 0.0
        
        for value in values:
            total += value
                
        if len(values) > 0:
            return [float(total / len(values))]
        else:
            return []
    
    def getValues(self,chrm, start, stop):
        return self.meanValue(self.treatment.getValuesInRange(chrm,start,stop))

#####

class Islands(DataBehaviour):
    
    def __init__(self,filename):
        self.treatment = BedFile(filename)
        self.getValues = self.treatment.getValuesInRange
    
    # bed islands uses avg(sum(  0/1([...]) , 0/1([...]), ...   )
    # does not normalise by bp
    def valuesBehaviour(self,chrm, values, start, end):
        if len(values)==0:
            return [0]
        else:
            return [1]
    

#####

from datastructures.AggregateTree import AggregateTree


#Tests for sufficient data in the row
def cpgMethValidRow(row):
    
    #Parameter for the minimum % of empty slices in a row (slices with no CpGs)
    cuttoff = 1.00
    counter = 0.00
    
    for i in row: 
        if i == None:
            counter +=1
    
    # Tests for a % of null greater than the cuttoff (i.e. if more than the cuttoff % of slices are empty return false)
    if counter / float(len(row)-1) >= cuttoff:
        return False
    else:
        return True


class CpGMethPercent(DataBehaviour):
    
    #Takes the CpG format input, plus left/right
    def __init__(self,filename,sample=0):
        assert sample.lower() in ["l","r","left","right", "absdiff", "reldiff", "folddiff"]
        
        ### treatment is THE TREE FOR CpGs - allows rapid searching.
        self.treatment = AggregateTree(filename)
        ### What is sample lower?
        self.sample = sample.lower()
        self.valuesBehaviour = missingValuesDontCount
    
    # For a window?
    def methPercentage(self,values):
        lmethTotal = 0
        lunmethTotal = 0
        rmethTotal = 0
        runmethTotal = 0
  
        for (lmeth,lunmeth,rmeth,runmeth) in values:
            lmethTotal += lmeth
            lunmethTotal += lunmeth
            rmethTotal += rmeth
            runmethTotal += runmeth
            
        lmeth = float(lmethTotal)
        lunmeth = float(lunmethTotal)
        
        # arbitrary cutoff for regions with very few reads
        if lmeth + lunmeth <= 10.0:
            return []
        
        rmeth = float(rmethTotal)
        runmeth = float(runmethTotal)
        
        # arbitrary cutoff for regions with very few reads
        if rmeth + runmeth <= 10.0:
            return []
        
        #Adds an ~0 value to each total (removes a division by 0 error)
        lmeth += 0.0000000001
        lunmeth += 0.0000000001
        rmeth += 0.0000000001
        runmeth += 0.0000000001
        
        # Calculates the percentages and differences
        lpercentage = lmeth / (lmeth+lunmeth)
        rpercentage = rmeth / (rmeth+runmeth)
        absDiff = rpercentage - lpercentage
        relDiff = absDiff/lpercentage
        foldDiff = math.log(rpercentage,2) - math.log(lpercentage,2)
        
        if relDiff > 5:
            relDiff = 5
        elif relDiff < -5:
            relDiff = -5
        if foldDiff > 5:
            foldDiff = 5
        elif foldDiff < -5:
            foldDiff = -5
            
        if self.sample.startswith("l"):
            return [lpercentage]
        elif self.sample.startswith("rel"):
            return [relDiff]
        elif self.sample.startswith("a"):
            return [absDiff]
        elif self.sample.startswith("f"):
            return [absDiff]
        else:
            return [rpercentage]
    
    # get the cpg meth percentage for a given range (chrm,start,stop)   
    def getValues(self,chrm, start, stop):
        return self.methPercentage(self.treatment.getValuesInRange(chrm,start,stop))
    
    # Tests for an appropriate R data range 
    def heatmapHasNegativeValues(self):
        
        if self.sample=="absdiff":
            return True
        elif self.sample=="reldiff":
            return True
        else:
            return False
    
    def heatmapUpperLowerBounds(self):
        if self.sample=="absdiff":
            return -1,1
        elif self.sample=="reldiff":
            return -5,5
        else:
            return None, None

class CpGMethPercentDifference(DataBehaviour):
    def __init__(self,filename):
        self.treatment = AggregateTree(filename)
        self.valuesBehaviour = missingValuesDontCount
    
    def methPercentage(self,values):
        lmethTotal = 0
        lunmethTotal = 0
        rmethTotal = 0
        runmethTotal = 0
  
        for (lmeth,lunmeth,rmeth,runmeth) in values:
            lmethTotal += lmeth
            lunmethTotal += lunmeth
            rmethTotal += rmeth
            runmethTotal += runmeth
            
        
        lmeth = float(lmethTotal)
        lunmeth = float(lunmethTotal)
        
        # arbitrary cutoff for regions with very few reads
        if lmeth + lunmeth <= 10.0:
            return []
        
        rmeth = float(rmethTotal)
        runmeth = float(runmethTotal)
        
        # arbitrary cutoff for regions with very few reads
        if rmeth + runmeth <= 10.0:
            return []
        
        lpercentage = lmeth / (lmeth+lunmeth)
        rpercentage = rmeth / (rmeth+runmeth)
        
        return [rpercentage - lpercentage]
    
    # get the cpg meth percentage for a given range (chrm,start,stop)   
    def getValues(self,chrm, start, stop):
        return self.methPercentage(self.treatment.getValuesInRange(chrm,start,stop))


#####

from sequence.genome import MemoryGenome as Genome

class GCContent(DataBehaviour):
    def __init__(self,build):
        self.genome = Genome(genomeBuild = build)
        
        self.valuesBehaviour = missingValuesDontCount
        self.chromosomeEnds = ChromosomeEnds(build)
    
    def gcPercent(self,sequence):
        Ns = sequence.count('N')
        
        Gs = sequence.count('G')
        Cs = sequence.count('C')
        
        effectiveSeqLen = len(sequence) - Ns
        
        assert effectiveSeqLen >= 0
        
        if effectiveSeqLen == 0:
            return None # not a valid region to work out the gc percentage of (all N's)
        else:
            # count the GC content in the sequence
            return float(Gs+Cs)/float(effectiveSeqLen)
        
    # calculates gc% content of the region
    def getValues(self,chrm, start, stop):
        if start > self.chromosomeEnds[chrm] or stop <= 0:
            return [] # not a valid region to work out the percentage of (entirely off chromosome)
        # retrieve the sequence for start -> end
        return [self.gcPercent(self.genome.getSequence(chrm,(max(start,0),min(stop,self.chromosomeEnds[chrm]))).upper())]      

#####

class CpGContent(DataBehaviour):
    def __init__(self,build):
        self.genome = Genome(genomeBuild = build)
        
        self.valuesBehaviour = missingValuesDontCount
        self.chromosomeEnds = ChromosomeEnds(build)
    
    # calculates cpg content of the region
    def cpgPercent(self,sequence):
        
        assert len(sequence) > 0
        
        Ns = sequence.upper().count('N')
        
        cpgs = sequence.upper().count('CG')
        
        effectiveSeqLen = len(sequence) - Ns
        
        if effectiveSeqLen == 0: # all N's
            return None
        else:
            cpgContent = float(cpgs) / float(effectiveSeqLen)
    
        return cpgContent
        
    def getValues(self,chrm, start, stop):
        if start > self.chromosomeEnds[chrm] or stop <= 0:
            return [] # not a valid region to work out the percentage of (entirely off chromosome)                
        # retrieve the sequence for start -> end
        return [self.cpgPercent(self.genome.getSequence(chrm,(max(start,0),min(stop,self.chromosomeEnds[chrm]))).upper())]     

#####

class CpGRatio(DataBehaviour):
    def __init__(self,build):
        self.genome = Genome(genomeBuild = build)
        
        self.valuesBehaviour = missingValuesDontCount
        self.chromosomeEnds = ChromosomeEnds(build)
    
    def cpgRatioCalc(self,chrm,start,stop,sequence):
        
        assert len(sequence) >= 0, sequence + ", " + chrm + ":" + str(start) +"-" + str(stop)
        
        Ns = sequence.upper().count('N')
        Gs = sequence.upper().count('G')
        Cs = sequence.upper().count('C')
        
        cpgs = sequence.upper().count('CG')
        
        effectiveSeqLen = len(sequence) - Ns
        
        if effectiveSeqLen == 0: # all N's
            return None
        else:                    
            if cpgs>0:
                cpgRatio = (float(cpgs) * effectiveSeqLen) / float(Gs*Cs) # from Schubeler (http://www.nature.com/ng/journal/v39/n4/full/ng1990.html)
            else:
                cpgRatio = 0.0 # avoid div by 0 error if Gs+Cs == 0, shortcut cpgs == 0 then cpgRatio = 0
    
        return cpgRatio
        
    def getValues(self,chrm, start, stop):
        if start > self.chromosomeEnds[chrm] or stop <= 0:
            return [] # not a valid region to work out the percentage of (entirely off chromosome)
        
        # retrieve the sequence for start -> end
        ratio = self.cpgRatioCalc(chrm,start,stop,self.genome.getSequence(chrm,(max(start,0),min(stop,self.chromosomeEnds[chrm]))).upper())
        
        return [] if ratio == None else [ratio]

#####

from datastructures.cintervaltree import Interval

#assumes all intervals have same chromosome
def mergeInts(intervals):      
        #Sorts the list
        intervals.sort(key=lambda x: (x.start, x.end))
       
        mergedIntervals = []
         
        #merges the list
        previous = None
        for current in intervals:
           
            if previous == None:
                previous = current
                continue

            # adjacent intervals should be merged, not just overlapping
            elif (current.end >= previous.start and current.start <= previous.end):
                start = min(current.start, previous.start)
                end = max(current.end, previous.end)
               
                previous = Interval(start,end)
               
            else:
                mergedIntervals.append(previous)
                previous = current
       
        mergedIntervals.append(previous)
       
        return mergedIntervals



class Coverage(DataBehaviour):
    
    def __init__(self,filename,build):
        self.treatment = BedFile(filename)
        self.getValues = self.treatment.getIntervalsInRange
        self.chromosomeEnds = ChromosomeEnds(build)
    
    def mergeIntervals(self,intervals):      
        return mergeInts(intervals)

    # percentage of coverage (by merging intervals and then seeing what percentage of BP is covered)
    def valuesBehaviour(self,chrm, intervals, start, end):
        if start > self.chromosomeEnds[chrm] or end <= 0:
            return [] # not a valid region to work out the percentage of
        
        # at least some of the slice must overlap
        if start < 0:
            start = 0 # trim if overlapping on left
        
        if end > self.chromosomeEnds[chrm]:
            # trim if overlapping on right
            end = self.chromosomeEnds[chrm]
        
        # rare but can happen if off end of chromosome
        if end - start == 0:
            return []
        
        self.mergeIntervals(intervals)
        bps = 0
        for interval in intervals:
            # intervals might extend outside this window, only count the portion that is internal to this window
            s = max(interval.start, start)
            e = min(interval.end,  end)
            bps += e-s
        
        assert bps <= end-start,  "More than 100% base usage"
        return [100.0 * float(bps) / float(end-start)]
    

