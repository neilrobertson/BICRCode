'''
Created on 25 Jul 2012

@author: mcbryan
'''

import collections
from itertools import izip
import csv
import math
        
# composite gene class
class CompositeGene():
    keys = None # keys
    
    def __init__(self, name, regionBehaviour,dataBehaviour,avgFunction,validRow):
        
        # Stores parameters
        # name = name of file
        # regionBehaviour = object for gene names, ids, positions etc
        # dataBehaviour = the data tree
        self.name = name # name of composite gene
        self.regionBehaviour = regionBehaviour
        self.dataBehaviour = dataBehaviour
        self.avgFunction = avgFunction
        self.validRow = validRow
        
        #Creates, dictionaries to store the slices, index is slice, list entry is each feature.
        #The second dictionary stores the lengths of each feature at each slice (as they differ by FEATURE. Used for normalising) 
        self.composite = collections.defaultdict(list) # dictionary of lists, key = index in slices
        self.slicewidths = collections.defaultdict(list) # as above
        
        self.fulloutput = [] # this is where the entire composite gene details are stored
        
        self.fulloutputsaved = None
    
    # for updating our composite structures with new values
    # Takes a feature (e.g. a gene) and slices, then updates dictionaries
    def addToComposite(self, obj):
        
        # gets the object containing gene names, ids, etc from Tonys parsed list, for the correct chromosome
        chrm = self.regionBehaviour.getChr(obj)
        
        # slices the gene:
        # Calls slicer (in regions), ultimately in Baselib - gene slicer.
        # Returns parallel lists of keys and slices, each entry is a string. REALLY EACH IS IN A SLICES OBJECT.
        keys,slices = self.regionBehaviour.slicer(obj)
        
        # only updates the keys if they look more complete
        ## Gets the keys list from the constructor, ~checks that all of the keys in the slice keys are the same as those in the constructor slice keys (initially these are empty)
        if CompositeGene.keys == None or len(CompositeGene.keys)<len(keys):
            CompositeGene.keys = keys
        
        # note that this is done in key order
        assert len(keys)==len(slices), "Make sure the list of keys is the same length as list of slices"
        
        rowvalues = {}
        
        # start full output with the id of the row (gene)
        ##This is simple a string of the feature (e.g. gene name)
        row = [str(obj)]
        #Gets the parallel key and slice information using izip (for each slice)
        for key,slice in izip(keys,slices):
            assert key == slice.sliceid, "Make sure the keys and sliceids are in the same order"
            
            # get raw data in the range slice.start -> slice.end
            ## Uses the data tree to pull raw data based on the slice start and end.
            values = self.dataBehaviour.getValues(chrm, slice.start, slice.end)
            
            # calculate value(s) for the slice
            ## Gets the values:
            ## Uses ultimately the L/R thing (this occured downstream.
            slicevalue = self.dataBehaviour.valuesBehaviour(chrm, values, slice.start, slice.end)
            
            rowvalues[slice.sliceid]=(slicevalue,slice.end-slice.start)
            
            # add a row for full output
            ## Row was the gene name earlier, this is the output!
            ## Gets the average... using the average function object, from right at the start.
            row.append(self.avgFunction(slicevalue))
            
        assert len(rowvalues)==len(slices),str(keys) + str(slices)
        
        if self.validRow(row):
            # Appends the slice VALUES (i.e. the final data, for each feature) to the full output.
            # This is in the format: feature name, slice 1 value, slice 2 value, slice 3 value, etc
            self.fulloutput.append(row)
            
            for sliceid in rowvalues:
                # add to composite gene
                ## Updates: for the current slice, gets the list from composite, then appends the values to the list (i.e. a new gene)
                ## Does the same for the lengths....
                slicevalue,slicelength = rowvalues[sliceid]
                self.composite[sliceid].extend(slicevalue)
                self.slicewidths[sliceid].append(slicelength)
        
            # this should always be the case after we have added the slices
            assert len(slices) <= len(self.composite),"Object: "+str(obj)+", Len(slices) ["+str(len(slices))+"] != len(self.composite) [" + str(len(self.composite))+"]"
    
    #Iterates for each vailidated feature (e.g. genes, exons etc)
    def buildComposite(self, regionsToUse):
        for obj in regionsToUse:
            # add the new values to the composite for those slices
            #Passes the feature, to be sliced and added to the dictionaries of lists
            self.addToComposite(obj)
    
    
    def writeFulloutput(self,filename):
        with open(filename,"w") as f:
            csvfile = csv.writer(f,delimiter="\t")
            
            header = [self.name]
            header.extend(self.keys)
            csvfile.writerow(header)
            
            csvfile.writerows(self.fulloutput)
            
            self.fulloutputsaved = filename
            
    def heatmapCommand(self, dataSource):
        slicewidths = [0.0] # prepopulate with 0 as breaks will be number of widths + 1
        for k in CompositeGene.keys:
            slicewidths.append(self.avgFunction(self.slicewidths[k]))
        
        scalingfactor = (len(slicewidths) - 1) / math.fsum(slicewidths) 
        
        breaks = []
        
        runningtotal = 0
        for width in slicewidths:
            runningtotal += width
            breaks.append(runningtotal * scalingfactor)
        
        commands = []
        
        commands.append("set.seed(0)")
        commands.append("x <- read.table (file='"+self.fulloutputsaved+"', header=TRUE, sep='\\t', quote='\\'', dec='.', fill=FALSE, comment.char='#', row.names=1,  na.strings = 'NA', nrows = -1, skip = 0, check.names = TRUE, strip.white = FALSE, blank.lines.skip = TRUE)")
        commands.append("library(gplots)")
        commands.append("pdf(file='"+self.fulloutputsaved+".pdf',width=4.5,height=9)")
            
        #Tests for the appropriate R colour scale
        if dataSource.heatmapHasNegativeValues():
            a,b = dataSource.heatmapUpperLowerBounds()  
            commands.append("heatmap.2(data.matrix(x),key=FALSE, Colv=FALSE, Rowv=FALSE, dendrogram='none',cexCol = 0.25, col=colorpanel(n=40,low='blue',mid = 'white', high='red'), breaks=seq(" + str(a) + "," + str(b) + ", length = 41), trace='none', labRow=' ', sepwidth=c(0,0), colWidths = c("+",".join(str(i) for i in breaks)+"))")  
        else:
            commands.append("heatmap.2(data.matrix(x),key=FALSE, Colv=FALSE, Rowv=FALSE, dendrogram='none',cexCol = 0.25, col=colorpanel(n=40,low='white',high='red'), breaks=sort(kmeans(Filter(function(x){x>=0}, c(data.matrix(x))),iter.max=50,centers=41)$centers), trace='none', labRow=' ', sepwidth=c(0,0), colWidths = c("+",".join(str(i) for i in breaks)+"))")  
        
        commands.append("dev.off()")
        
        return "\n".join(commands)
    