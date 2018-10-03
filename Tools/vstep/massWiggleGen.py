#!/usr/bin/env python

import os
import sys
import getopt
import csv
import subprocess
import Queue
import threading
from datastructures.cintervaltree import *
import pylab
import numpy
import math

from genemapping import Ensembl
import genemapping.geneslicer
from csvfile.genelist import GeneList

def makeDirectory(dir):
    if not os.access(dir, os.R_OK):
        print "making directory", dir
        os.mkdir(dir)
        return dir

if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:v:e:o:p:f:w:", [])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    # add executing directory as part of path
    sys.path.append(sys.path[0])
    
    vstepWidth = 200 #default
    numbSlices = 20 # default

    for o, a in opts:
        # need a vstep width
        if o=="-s":
            vstepWidth = int(a)
        # need a vstep file
        elif o == "-v":
            vstepFile = a
        # need a list of genes -> coords
        elif o == "-e":
            genesFileName = a # eg : "genes-and-exons-human-NCBI36.csv"
            print "Loading gene mapping"
            genesmapping = Ensembl.GenesMapping(genesFileName)
        # need an output folder
        elif o == "-o":
            outputFolder = a
        # need a number of slices
        elif o == "-p":
            numbSlices = a
        elif o == "-f":
            friendlyGenesNames = Ensembl.FriendlyGeneNames(a)
        elif o == "-w":
            webroot = a
        
    # need lists of genes
    geneListLocations = args

    # create output folder
    makeDirectory(outputFolder)

    mountFolder = os.path.expanduser("~/mount")

    stdout_lock = threading.Lock()
    figure_lock = threading.Lock()
    analysis_list_lock = threading.Lock()

    analysedgenes = {} # map of set -> [] of genes

    # a bit of parallelisation
    def dowork(outputFolder, geneListName, gene, chr, strand, genestart, genestop):
        
        stdout_lock.acquire()
        print geneListName, " - Working on : ",  gene, chr, strand, genestart, genestop
        stdout_lock.release()

        # 5kbp either side of the gene
        extendedstart = genestart - 5000
        extendedstop  = genestop + 5000
        
        wiggleLocation = outputFolder + "/" + geneListName + "/" + gene + ".wig"
        
        # vstepToWiggle
        #./vstep2wiggle.sh vstepFile vstepWidth chr start stop
        wigglegenerator = subprocess.Popen(args=[mountFolder+"/repository/shared/baselib/vstep2wiggle.sh", vstepFile, str(vstepWidth), chr, str(extendedstart), str(extendedstop)],  stdout=open(wiggleLocation, "w"))
        wigglegenerator.wait()
        del wigglegenerator
        
        # read in the wiggle and draw some charts with it
        
        # build small interval tree from the wiggle file
        tree =  IntervalTree()
        wiggleFile = open(wiggleLocation, "r")
        wiggle = csv.reader(wiggleFile, delimiter="\t")
        wiggle.next() # skip first line (header)
        for row in wiggle:
            (chr, start, stop, value) = row
            start = int(start)
            stop = int(stop)
            value = int(value)
            tree.insert_interval(Interval(start, start+vstepWidth, value=value))
        wiggleFile.close()
        del wiggleFile
        del wiggle
        
        # slice gene
        keys,  slices = geneslicer.sliceGene((chr, strand, genestart, genestop), numbSlices, 5000, 1000, 5000, 1000)
        
        # obtain values along the gene
        values = {}
        for slice in slices:
            intervals = tree.find(min(slice.start, slice.end), max(slice.start, slice.end))
            
            alignments = []
            for interval in intervals:
                alignments.append(interval.value)
    
            values[slice.sliceid] = 0 if len(alignments)==0 else math.fsum(alignments) / float(len(alignments))
            del intervals
            del alignments
        del tree
        
        # output plot data
        plotData = outputFolder + "/" + geneListName + "/" + gene + ".plot"
        plotdataFile = open(plotData, "w")
        plotdata = csv.writer(plotdataFile, delimiter="\t")
        for key in keys:
            plotdata.writerow([key, values[key]])
        plotdataFile.close()
        del plotdataFile
        del plotdata
        
        # plot figure
        plotLocation = outputFolder + "/" + geneListName + "/" + gene + ".png"
        figure_lock.acquire()
        pylab.figure(1)
        pylab.plot([values[key] for key in keys],'.-',label=gene )
        pylab.xticks(range(len(keys)), keys,  rotation=45)
        pylab.xlabel(gene)
        pylab.ylabel('Avg Expression')
        pylab.ylim(ymin=0)
        pylab.savefig(plotLocation,  dpi=640/8)
        pylab.clf()
        figure_lock.release()
        
        analysis_list_lock.acquire()
        if not geneListName in analysedgenes:
            analysedgenes[geneListName] = {}
        analysedgenes[geneListName][gene] = (gene, chr, strand, genestart, genestop, wiggleLocation, plotLocation,  plotData)
        analysis_list_lock.release()
        
        del values
        del keys

    # set up worker threads
    cores = 15
    wigglegenerators = Queue.Queue(cores)

    def worker():
        while True:
            item = wigglegenerators.get()
            dowork(*item)
            wigglegenerators.task_done()
            del item
            
    for i in range(cores):
        t = threading.Thread(target=worker)
        t.setDaemon(True)
        t.start()

    # load gene lists
    geneLists = {}
    for geneListLocation in geneListLocations:
        genes = GeneList(geneListLocation)
        geneListName = genes.getFriendlyName()
        makeDirectory(outputFolder+"/"+geneListName)
        print "Loading gene set:", geneListName

        geneLists[geneListName] = genes

    # for each gene list do some wiggling / charting
    for geneList in geneLists:
        genes = geneLists[geneList]
        # for each gene
        for gene in genes:
            # look up coordinates of gene
            if gene in genesmapping:
                chr, strand, start, stop = genesmapping[gene]
                
                wigglegenerators.put((outputFolder, geneList, gene, chr, strand, start, stop))
    
    # wait till its done
    wigglegenerators.join()
    
    # make an index file for each gene list
    for geneList in geneLists:
        htmlFile = open(outputFolder + "/" + geneList + "/" + "index.html", "w")
        
        print >>htmlFile,  "<html><body>"
        
        print >>htmlFile,  "<table border=\"1\">"
        
        for individualgene in geneLists[geneList]:
            if individualgene in genesmapping:
                (gene, chr, strand, genestart, genestop, wiggleLocation, plotLocation, plotData) = analysedgenes[geneList][individualgene]
                
                wiggleLocation = wiggleLocation[wiggleLocation.rfind("/")+1:]
                plotLocation = plotLocation[plotLocation.rfind("/")+1:]
                plotData = plotData[plotData.rfind("/")+1:]
                
                # http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=chr4:114661572-114867995&hgt.customText=http://130.209.136.101//pythonapps/testbed/files/2375_46.wig.txt
                
                print >>htmlFile,  "<tr>"
                print >>htmlFile, "<td>"+friendlyGenesNames[gene]+"</td>"
                print >>htmlFile,  "<td><b><a href=\"http://may2009.archive.ensembl.org/Homo_sapiens/Gene/Summary?g="+gene+"\">"+gene+"</a></b></td>"
                print >>htmlFile,  "<td>"+chr+":"+str(genestart)+"-"+str(genestop)+", "+strand+" strand</td>"
                print >>htmlFile,  "<td><p><a href=\""+wiggleLocation+'">Wiggle</a></b></p>\
                <p><a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=' \
                                                                                                                +chr+":"+str(genestart-5000)+"-"+str(genestop+5000) \
                                                                                                                +'&hgt.customText='+webroot+"/" +geneList +"/"+wiggleLocation+'">UCSC Browser</a></p>\
                </td>'
                print >>htmlFile,  "<td></td>"
                print >>htmlFile,  '<td><a href="' + plotLocation + '"<img src="'+plotLocation+'" height="200" alt="Plot"/></a> \
                                                      <a href="'+plotData+'">PlotData</a></b></td>'
                print >>htmlFile,  '<td></td>'
                print >>htmlFile,  "</tr>"
            
        print >>htmlFile,  "</table>"
            
        print >>htmlFile,  "</html></body>"
