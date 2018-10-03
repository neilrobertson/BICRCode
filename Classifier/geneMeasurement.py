#!/usr/bin/env python

import os
import sys

# add executing directory as path
sys.path.append(sys.path[0])
sys.path.append(os.path.expanduser("~/mount/repository/shared/baselib"))

import csv
import getopt
import numpy
import math

from vstep.treatment import Vstep as VstepTreatment
from bed.treatment import Bed as BedTreatment

from genemapping.pointlist import Pointlist

from genemapping.Ensembl import GenesMapping as EnsemblGenesMapping
from genemapping.Affy import GenesMapping as AffyGenesMapping

from genemapping.chrmEnds import ChromosomeEnds

from csvfile.genelist import GeneList as GenesToUse
import genemapping.geneslicer as geneslicer

def getGeneChr(gene):
    (chr, strand, start, end) = gene
    return chr
    
def getPointChr(point):
    (chr, point) = point
    return chr

def getGeneSlice(gene, padding):
    (chr, strand, start, end) = gene
    
    assert strand == '+' or strand == '-'
    if strand == '+':
        return (start-padding, end+padding)
    else:
        return (start-padding, end+padding)

def getTssPaddingGeneSlice(gene, padding):
    (chr, strand, start, end) = gene
    
    assert strand == '+' or strand == '-'
    if strand == '+':
        return (start-padding, end)
    else:
        return (start, end+padding)
        

def getTSSSlice(gene, padding):
    # completely ignore padding parameter
    (chr, strand, start, end) = gene
    
    uppadding = 5000
    downpadding = 1000
    
    assert strand == '+' or strand == '-'
    if strand == '+':
        return (start-uppadding, start+downpadding)
    else:
        return (start-downpadding, start+uppadding)
    
def getPointSlice(point, padding):
    (chr, point) = point
    return (point-padding, point+padding)

# treatment behaviours
# this is how it combines together observances seen on a single gene slice

def treatmentSum(chr, values, start, end):
    # sum = 0 if no values
    assert math.fsum(values)>=0.0,  "Sum negative"
    return [math.fsum(values)]

def treatmentSumByBP(chr, values, start, end):
    assert math.fsum(values)>=0.0,  "Sum negative"
    return [math.fsum(values) / float(end-start)]

# for island files
def treatmentYesNo(chr, values, start, end):
    if len(values)==0:
        return [0]
    else:
        return [1]

# doesnt affect the calculation if there are no results
def treatmentAvgNoZero(chr, values, start, end):
    assert math.fsum(values)>=0.0,  "Sum negative"
    if len(values)==0:
        return []
    else:
        return [math.fsum(values)/float(len(values))]

def mergeIntervals(intervals):
    for i in range(len(intervals)):
        for j in range(i, len(intervals)):
            if i == j:
                continue
            else:
                if (intervals[j].start < intervals[i].end and intervals[i].start < intervals[j].end):
                    intervals[i].start = min(intervals[i].start,  intervals[j].start)
                    intervals[i].end = max(intervals[i].end,  intervals[j].end)
                    intervals.remove(intervals[j])
                    mergeIntervals(intervals)
                    return

def treatmentPercentageOfBP(chr, intervals, start, end):
    if start > chromosomeEnds[chr] or end < 0:
        return [] # not a valid region to work out the percentage of
    
    # at least some of the slice must overlap
    if start < 0:
        start = 0 # trim if overlapping on left
    
    if end > chromosomeEnds[chr]:
        # trim if overlapping on right
        end = chromosomeEnds[chr]
    
    # rare but can happen
    if end - start == 0:
        return []
    
    mergeIntervals(intervals)
    bps = 0
    for interval in intervals:
        # intervals might extend outside this window, only count the portion that is internal to this window
        s = max(interval.start, start)
        e = min(interval.end,  end)
        bps += e-s
    
    assert bps <= end-start,  "More than 100% base usage"
    return [100.0 * float(bps) / float(end-start)]
    

# a composite gene class
class GeneSummary():
    surroundingDistance = 5000
    
    def __init__(self, name):
        self.name = name # name of composite gene
    
    def doSummary(self, genesToUse, treatment):
        #self.lengths = []
        for genename, gene  in genesToUse:
            # get gene id details
            chr = getChr(gene)
            
            start, end = getSlice(gene,  self.surroundingDistance)
            
            values = getIntervals(chr, start, end)
            
            result = treatmentBehaviour(chr, values, start, end)
            
            # the value for the gene would be here
            print genename ,"\t", result[0]

if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ds:v:a:b:e:t:i:p%:c:g:", [])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        print "Usage: main.py -s <vstepSize> -v/b/i/% <treatmentValuesFileName> -a/e/p <genesMappingFileName> [Space seperated list of gene id sets]"
        sys.exit(2)

    # add executing directory as part of path
    sys.path.append(sys.path[0])
    
    vstepWidth = 200 #default
    
    pointlist = False
    
    for o, a in opts:
        if (o=="-s"):
            vstepWidth = int(a)
        elif (o=="-d"):    
            debug = True
        elif (o=="-c"):
            chromosomeEnds = ChromosomeEnds(a)
        elif (o=="-g"):
            GeneSummary.surroundingDistance = int(a)

    for o, a in opts:
        if o == "-s":
            pass # we dealt with this already
        elif (o=="-d"):
            pass # we dealt with this already
        elif (o=="-c"):
            pass # we dealt with this already
        elif (o=="-g"):
            pass # we dealt with this already

        # gene mappings / point list settings (for how to deal with the arguments representing which genes / points to plot)
        elif o == "-a":
            genesFileName = a # eg : "HG-U133Plus2.csv"
            genesmapping = AffyGenesMapping(genesFileName)
            getChr = getGeneChr
            getSlice = getGeneSlice
        elif o == "-e":
            genesFileName = a # eg : "genes-and-exons-human-NCBI36.csv"
            genesmapping = EnsemblGenesMapping(genesFileName)
            getChr = getGeneChr
            getSlice = getGeneSlice
        elif o == "-t":
            # transcription start site
            genesFileName = a # eg : "genes-and-exons-human-NCBI36.csv"
            genesmapping = EnsemblGenesMapping(genesFileName)
            getChr = getGeneChr
            getSlice = getTSSSlice
        elif o == "-p":
            surroundingDistance = 100000
            pointlist = True
            getChr = getPointChr
            getSlice = getPointSlice

    
    
        # treatments (for what data we will be summarising)
        elif o == "-b":
            # bed uses avg(  sum([...]) , sum([...]), ...   )
            treatmentFileName = a
            treatmentBehaviour = treatmentSumByBP
            treatment = BedTreatment(treatmentFileName)
            getIntervals = treatment.getValuesOfOverlappingIntervals
        elif o == "-i":
            # bed islands uses avg(sum(  0/1([...]) , 0/1([...]), ...   )
            # does not normalise by bp
            treatmentFileName = a
            treatmentBehaviour = treatmentYesNo
            treatment = BedTreatment(treatmentFileName)
            getIntervals = treatment.getValuesOfOverlappingIntervals
        elif o == "-%":
            treatmentFileName = a
            treatmentBehaviour = treatmentPercentageOfBP
            treatment = BedTreatment(treatmentFileName)
            getIntervals = treatment.getOverlappingIntervals
        else:
            assert False,  "Unhandled option: " + o
    
    for genesToUseFileName in args:
        
        print "***" + genesToUseFileName + "***"
        geneSummary= GeneSummary(genesToUseFileName)

        if pointlist == True:
            genes = Pointlist(genesToUseFileName)
        else:
            genes = set()
            genesToUse = GenesToUse(genesToUseFileName)
            
            # build set of genes
            for gene in genesToUse:
                if not gene in genesmapping:
                    # Gene not found in coords mapping
                    continue
                else:
                    genes.add((gene, genesmapping[gene]))
        
        geneSummary.doSummary(genes, treatment)
