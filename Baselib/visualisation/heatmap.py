'''
Created on 23 Aug 2010

@author: mcbryan
'''
from datastructures.memoized import memoized
from datastructures.orderedset import OrderedSet
import numpy
import os.path
import csv
import sys
import rpy2.rpy_classic as rpy
#from rpy2.rpy_classic import *
import rpy2.robjects as robjects
import getopt
import collections
import itertools
import math
from rpy2.robjects.numpy2ri import numpy2ri
from csvfile.genelist import GeneList

rpy.set_default_mode(rpy.NO_CONVERSION)

na_warning_shown = False

def naBehaviour(arr,na_to_zeroes,remove_na):
    # can't have both at once
    assert not (na_to_zeroes and remove_na)
    if na_to_zeroes:
        return map(lambda x: "0.0" if x=="NA" else x, arr)
    if remove_na:
        for a in arr:
            if a == "NA":
                if not na_warning_shown:
                    print "I seen NA's - I'm removing these rows entirely.  Hopefully this is what you want."
                    na_warning_shown = True
                raise Exception()
    return arr


######
# Load the raw expressions for heatmaps etc
###### 
@memoized
def loadExpressionData(expressionsfile, na_to_zeroes = False, remove_na = True):
    print "Loading Expression Data...",
    with open(expressionsfile) as infile:
        expressions = csv.reader(infile, delimiter="\t")

        firstColSeen = False

        ids = []
        col_names = []
        exprs = []

        for row in expressions:
            if not firstColSeen:
                row[0]="Name"
                col_names = row[1:]
                firstColSeen = True
            else:
                # put rma data into array
                try:
                    exprs.append(map(float, naBehaviour(row[1:],na_to_zeroes,remove_na)))
                    ids.append(row[0])
                except Exception:
                    continue
        

        assert len(exprs)==len(ids)
        print "Done. Loaded "+str(len(ids))+" rows."
        return ids, col_names, exprs
    

def drawHeatMap(exprsFileLoc,outputfileloc,rowNameFunction,arrayColourFunction,array_phenotypes,genes,genesOnceEach=False,allGenes=False,rowScaling=False,onlyConditions=None,dendrogram="row",heatmapColours=("red","black","green"),allowSkipping=True,log=False,indataorder=False):
    
    ids,col_names,exprs = loadExpressionData(exprsFileLoc)
    
    if allGenes == True:
        set_of_ids = ids
    else:
        set_of_ids = genes
    
    # reorder the list of genes if we want it to be in order of the raw data
    if indataorder:
        data = list()
        for id in ids:
            if id in set_of_ids:
                data.append(id)                
        set_of_ids = data
    
    rowNumbers = collections.defaultdict(list)
    for i in range(0,len(ids)):
        rowNumbers[ids[i]].append(i)
    
    row_names_subset = []
    newexprs = []
    
    for id in OrderedSet(set_of_ids):
        for rownumber in rowNumbers[id]:
            row_names_subset.append(id)
            if onlyConditions != None:
                # only add the columns we need here
                expression = []                    
                for j in range(0, len(col_names)):
                    if array_phenotypes[col_names[j]] in onlyConditions:
                            expression.append(exprs[rownumber][j])
            else:
                expression = exprs[rownumber]
            newexprs.append(expression)
            
            # exit loop after adding first one if once each is on
            if genesOnceEach:
                break
    
    # column names updated here
    if onlyConditions != None:
        new_column_names = []
        for i in range(0, len(col_names)):
            if array_phenotypes[col_names[i]] in onlyConditions:
                new_column_names.append(col_names[i])
    else:
        new_column_names = col_names
        
    exprs_as_array= numpy.array(newexprs, float)
    
    if log:
        # small offset to prevent log2(0) errors
        exprs_as_array = numpy.log2(exprs_as_array+0.01)
    
    array_colours = map(arrayColourFunction, new_column_names)
    
    if len(exprs_as_array)<2:
        print "Not plotting Heatmap "+outputfileloc+" of only "+str(len(exprs_as_array))+" rows"
        return
    
    if len(exprs_as_array)>60000:
        # sanity check
        print "Not plotting Heatmap "+outputfileloc+" of too many ("+str(len(exprs_as_array))+") rows"
        return
    
    print "Plotting Heatmap "+outputfileloc+" with " + str(len(exprs_as_array))+" rows",
    
    if allowSkipping and os.path.isfile(outputfileloc):
        print "...skipping (already exists)"
        return
    else:
        print # just to make it line up properly
    

    # calculate page size
    numbGenes = len(exprs_as_array)
    genesPerPage = 100 # 100 genes per "page"
    pagesToOutputOn = -(-numbGenes / genesPerPage) # rounds up always

    # convert exprs to R format
    exprs_as_array = numpy2ri(exprs_as_array)
        
    

    rpy.r.library("fastcluster")
    # this is converted to python code from the heatmap.2.R script in Gplots
    # we just want to do it here so it clusters better
    if rowScaling == True:
        means = rpy.r.rowMeans(exprs_as_array,True)
        exprs_as_array = rpy.r.sweep(exprs_as_array, 1, means)
        sx = rpy.r.apply(exprs_as_array, 1, rpy.r.sd, True)
        exprs_as_array = rpy.r.sweep(exprs_as_array, 1, sx, "/")


    
    # get a list of complete cases (genes that we will keep)
    complete_cases = rpy.r.complete_cases(exprs_as_array).sexp
    
    # filter the row names by this list    
    row_names_subset = [d for d, s in itertools.izip(row_names_subset, complete_cases) if s]
    
    # remove na's which can't be plotted (using same criteria as complete cases
    robjects.r('removeNAs <- function(x) x[complete.cases(x),]')
    exprs_as_array = rpy.r.removeNAs(exprs_as_array)

    # heatmap colours
    rpy.r.library("gplots")  
    # heatmap colours
    assert len(heatmapColours)==2 or len(heatmapColours)==3
    if len(heatmapColours) == 2:
        a,b = heatmapColours
        heatmapColours = rpy.r.colorpanel(75,a,b)
    else:
        a,b,c = heatmapColours
        heatmapColours = rpy.r.colorpanel(75,a,b,c)
    
    
        
    # This is important
    # heatmap_2 uses the isTrue() function to determine if values are TRUE or FALSE rather than using ==
    # this means that it must be the actual TRUE logical, not just something that evaluates to TRUE
    # unfortunately rpy2 seems to have r.TRUE and r.FALSE set to 1 and 0 respectively which eval correctly
    # but aren't really the right types
    
    # R does have T and F objects which are 
    true = robjects.r['T']
    false = robjects.r['F']    
    


    rpy.r.pdf(outputfileloc, paper="special",
            pointsize=12,
            #width=29.7/2.54,
            width=60/2.54,
            height=20.9/2.54 * pagesToOutputOn,
            pagecentre=True)
    
    matrix = rpy.r.matrix([5, 4, 0, 1, 3, 2], byrow=True,  nrow=3)
    
    def maxLabelLength(labels):
        return max([len(label) for label in labels])
    
    marginfactor = 12.0
    
    rowlabels = map(rowNameFunction, row_names_subset)
    
    rowmargin = min(40,max(10, math.ceil(marginfactor*maxLabelLength(rowlabels)/math.sqrt(len(rowlabels)))))
    colmargin = 15
    
    rpy.r.heatmap_2(exprs_as_array,
                labRow=rowlabels,
                labCol=new_column_names,
                Colv=false if dendrogram == "row" or dendrogram == "none" else true,
                Rowv=false if dendrogram == "col" or dendrogram == "none" else true,
                symbreaks=true if rowScaling else false,
                margins=[colmargin, rowmargin],
                ColSideColors=array_colours,
                dendrogram = dendrogram, 
                col=heatmapColours,
                key=true,
                symkey=false,
                #density="none",
                trace="none",
                #cexRow=0.5, 
                #cexCol=0.5,
                cexRow=1, 
                cexCol=1,
                lhei = [0.1, 0.05, float(pagesToOutputOn)/2.0],
                lwid = [1, 4],
                lmat = matrix
                )
    
    rpy.r.dev_off()
     
    # write the output - commented out currently as although it works R has all the columns shifted by 1 so it requires manual
    # fixing, also, pollutes file directory if not needed   
    #rpy.r.write_table(exprs_as_array,outputfileloc+".source",sep="\t",col = new_column_names, row = map(rowNameFunction, row_names_subset))


if __name__ == "__main__":

    try:
        opts, args = getopt.getopt(sys.argv[1:], "e:o:", ["data=","output=","colours=","rowscaling=","dendrogram=","unique","log","genelist=","orderbydata"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    exprs = None
    outputfileloc = None
    heatmapColours = ("green","black","red")
    rowScaling = True
    dendrogram = "row"
    unique = False
    log = False
    genelist = []
    indataorder = False

    for o, a in opts:
        if o in ["-e","--data"]:
            exprs = a
        elif o in ["-o","--output"]:
            outputfileloc = a
        elif o in ["--colours"]:
            heatmapColours = tuple(a.split(","))
            assert 2 <= len(heatmapColours) <= 3
        elif o in ["--rowscaling"]:
            rowScaling = a.lower() in ("yes", "true", "t", "1")
        elif o in ["--dendrogram"]:
            assert a in ["row","column","both","none"]
            dendrogram = a
        elif o in ["--unique"]:
            unique = True
        elif o in ["--log"]:
            log = True
        elif o in ["--orderbydata"]:
            assert dendrogram in ["column","none"]
            indataorder = True
        elif o in ["--genelist"]:
            genelist = GeneList(a)

    assert exprs != None
    assert outputfileloc != None

    if len(genelist) != 0:
        print genelist
        print len(genelist)

    drawHeatMap(exprs,
                outputfileloc,
                lambda x : x,
                lambda x : "white",
                lambda x : "a",
                genelist,
                allGenes = True if len(genelist) == 0 else False,
                rowScaling=rowScaling,
                onlyConditions=None,
                dendrogram = dendrogram,
                heatmapColours = heatmapColours,
                genesOnceEach = unique,
                allowSkipping=False,
                log=log,
                indataorder=indataorder
                )

