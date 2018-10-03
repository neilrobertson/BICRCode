#!/usr/bin/env python

import sys

import os
import csv
import getopt
import numpy
import math

from Composite import CompositeGene

from Regions import AffyRegions, EnsemblRegions, PointRegions, BedRegions
print(sys.version)
from Data import Vstep,VstepIntermediateAverage,Bed,BedCount,BedValues,Coverage,Islands,CpGContent, GCContent, CpGRatio, CpGMethPercent, CpGMethPercentDifference, BedWithControl, cpgMethValidRow

def getAvg(values):
    if len(values)==0:
        return None
    return math.fsum(values) / float(len(values))

def getStdDev(windowValues):
    return numpy.std(windowValues)

def getSEM(windowValues):
    return getStdDev(windowValues) / numpy.sqrt(len(windowValues))

def validRow(row):
    return True

if __name__ == "__main__":
    
    # c - assembly
    
    # treatment to plot details
    # b - bed <bed file>
    # i - islands <islands bed file>
    # v - vstep <vstep file> # s - vstep width <integer>
    # % - percentage <percentage coverage bed file>
    
    # feature list modifiers

    # a - affy <affy mapping file>
    # e - ensembl
    # p - points
    
    # extend - used for bed files only
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", [
        "assembly=", 
        "annotation=",
        
        # what are we plotting
        "bed=", "bedCntrl=", "bedCount=","bedValues=",
        "vstep=", "vstep-zeroes=",  "vstep-intermediate=", "vstep-zeroes-intermediate=",  "vstepwidth=",
        "cpg-meth-percent-difference=",
        "cpg-meth-percent=","cpg-meth-percent-sample=",
        "islands=",
        "coverage=",
        
        # built in metrics (gc content, cpg content, cpg ratio
        "gccontent",
        "cpgcontent",
        "cpgratio",
        
        # what are our features (affy requires a mapping file as well)
        "affy=", 
        "ensembl","tsses", "exons=", "min-exons=",
        "bedList",
        "points",
        "fulloutput=",
        "summary=",
        
        # interval distances etc
        "ud=",
        "ui=",
        "dd=",
        "di=",
        "chunks=",
        "extend=",
        
        
        "printEntries"
        ])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        print "Usage: main.py  [Space seperated list of gene id sets]"
        sys.exit(2)
    
    vstepWidth = None #default
    
    exons = False
    minExons = 3
    numbExons = None
    fulloutput = None
    tsses = False
    build = None # default
    cpgMethPercentSample = None
    annotation = "ncbi36.1"
    
    printEntries = False 
    
    ud = None
    ui = None
    dd = None
    di = None
    chunks = None
    extend = None
    summary = sys.stdout
    
    # parameters that need to be parsed first
    for o, a in opts:
        if (o=="--vstepwidth"):
            vstepWidth = int(a)
        elif (o=="--assembly"):
            build = a
        elif (o=="--annotation"):
            annotation = a
        elif (o=="--exons"):
            exons = True
            numbExons = int(a)
        elif (o=="--min-exons"):
            minExons = int(a)
        elif (o=="--fulloutput"):
            fulloutput = a
        elif (o=="--tsses"):
            tsses = True
        elif (o=="--ud"):
            ud = int(a)
        elif (o=="--ui"):
            ui = int(a)
        elif (o=="--dd"):
            dd = int(a)
        elif (o=="--di"):
            di = int(a)
        elif (o=="--chunks"):
            chunks = int(a)
        elif (o=="--extend"):
            extend = int(a)
        elif (o=="--cpg-meth-percent-sample"):
            assert a.lower() in ["l","r","left","right", "absdiff", "reldiff", "folddiff"]
            cpgMethPercentSample = a.lower()
        elif (o=="--summary"):
            summary = open(a,"w")
        elif (o=="--printEntries"):
            printEntries = True

    # load in data
    for o, a in opts:        
        # -------------- Regions to be interrogated settings ---------------------

        # gene mappings / point list settings (for how to deal with the arguments representing which genes / points to plot)
        if o == "--affy":
            regionSettings = AffyRegions(a, chunks, ud, ui, dd, di)

        #### Loads in gene names, ids, postitions, intron, exons etc. This uses the upstream dimension (ud), upstream interval (ui)
        #### Doesnt load the input (e.g. genes.nochange) files yet.
        elif o == "--ensembl":
            assert build != None
            if build != "hg18" and annotation == "ncbi36.1":
                assert False, "Trying to use default ncbi36.1 with non hg18 genome"
            regionSettings = EnsemblRegions(build, annotation, exons, minExons, numbExons, tsses, chunks, ud, ui, dd, di)
        
        elif o == "--points":
            assert build != None
            regionSettings = PointRegions(build, ud, ui, dd, di)
        
        elif o == "--bedList":
            assert build != None
            regionSettings = BedRegions(build, chunks, ud, ui, dd, di)

        # -------------- Raw data sources  ---------------------
        
        elif o == "--vstep":
            dataSource = Vstep(a, vstepWidth)
            
        elif o == "--vstep-zeroes":            
            dataSource = Vstep(a, vstepWidth, zeroes=True)
            
        elif o == "--vstep-intermediate":            
            dataSource = VstepIntermediateAverage(a, vstepWidth, getAvg)

        elif o == "--vstep-zeroes-intermediate":
            dataSource = VstepIntermediateAverage(a, vstepWidth, getAvg, zeroes=True)

        elif o == "--bed":
            dataSource = Bed(a, extend=extend)
        
        elif o == "--bedCount":
            dataSource = BedCount(a, extend=extend)
    
        elif o == "--bedValues":
            dataSource = BedValues(a, extend=extend)  

        elif o == "--islands":
            dataSource = Islands(a)
        
        #### Inititalises the CpG tree, from the CpG aggregate file ####
        elif o == "--cpg-meth-percent":
            assert cpgMethPercentSample != None
            dataSource = CpGMethPercent(a,cpgMethPercentSample)
            validRow = cpgMethValidRow
        
        elif o == "--cpg-meth-percent-difference":
            dataSource = CpGMethPercentDifference(a)
        
        elif o == "--gccontent":
            dataSource = GCContent(build)

        elif o == "--cpgcontent":
            dataSource = CpGContent(build)
            
        elif o == "--cpgratio":
            dataSource = CpGRatio(build)

        elif o == "--coverage":
            dataSource = Coverage(a,build)
    
    # stuff that needs to be done after primary data is loaded
    for o,a in opts:
        if o == "--bedCntrl":
            dataSource = BedWithControl(a,dataSource,extend=extend)
    
    compositeGenes = []
    
    #Iterates for each file in the --ensembl option
    for genesToUseFileName in args:
        
        print >> sys.stderr, "***" + genesToUseFileName + "***"
        
        # Region settings is the object containing dictionaries of gene Ids, positions etc.
        # Data source is the CpG loci tree
        # getAvg passes the object, getAvg, which can later be used to get averages
        # Creates a new composite gene object, for the current file.
        compositeGene = CompositeGene(genesToUseFileName,regionSettings,dataSource,getAvg,validRow)
    
        #Passes the file to the gene list class (ultimately) - gets a list of genes in the file.
        items = regionSettings.regionIterator(genesToUseFileName)
    
        #Creates a new list, testing that the genes in the current input file are in the regions list.
        # build set of valid items
        validItems = list()
        for item in items:
            if regionSettings.valid(item):
                validItems.append(item)
        
        print >> sys.stderr, "Valid items: ",str(len(validItems))

        # Does the work....?
        # Calls build composite> , which iterates through the list
        # At each iteration calls addtocomposite which:
        # 
        compositeGene.buildComposite(validItems)
        
        # Tests that the user wants a full output...
        # Includes keys etc - you should test this....
        if fulloutput != None:
            if ".csv" in fulloutput:
                # we have a full qualified file name we should use
                compositeGene.writeFulloutput(fulloutput)
            else:
                # we have a pattern to use for naming the composite genes
                compositeGene.writeFulloutput(fulloutput+"."+compositeGene.name[compositeGene.name.rfind("/")+1:]+".csv")
        
        # Appends the files composite results to the GLOBAL (all files) table.
        compositeGenes.append(compositeGene)

    # print composite profiles for all input sets

    # header
    print >> summary, "\t",
    for gene in compositeGenes:
        print >> summary, gene.name[gene.name.rfind("/")+1:] + "\t",
    print >> summary, ""
    
    # details
    for k in CompositeGene.keys:
        print >> summary, k + "\t", 
        for gene in compositeGenes:
            print >> summary, str(getAvg(gene.composite[k])) + "\t",
            # Tests for negative values
        print >> summary,""
   
    if printEntries:
        for k in CompositeGene.keys:
            print getAvg(gene.slicewidths[k])
            print >> summary, "len "+ k + "\t", 
            for gene in compositeGenes:
                print >> summary, str(len(gene.composite[k])) + "\t",
            print >> summary,""  
    
    # command to plot the heatmap in R, variable determines two or three colour.
    if fulloutput!=None:
        for compositeGene in compositeGenes:
            print
            print compositeGene.heatmapCommand(dataSource)
