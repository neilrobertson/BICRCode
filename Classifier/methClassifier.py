#!/usr/bin/env python

import os
import sys
import csv
import getopt
import math
from genemapping import Ensembl
from bed.treatment import Bed, ExtendedBed
from csvfile.indexedcsv import IndexedCSV
from affy.NetAffxAnnotation import NetAffxAnnotation

try:
    opts, args = getopt.getopt(sys.argv[1:], "i:o:a:ms", [])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    sys.exit(2)

# takes in a csv file of (point) coordinates and tells us some stuff about them

affyComparisonFile = None
onlyMethProbesWithAffyProbes = False
allRows = True

for o, a in opts:
    if o=="-i":
        infile = a
    elif o=="-o":
        outfile = a
    elif o=="-a":
        affyComparisonFile = a
    elif o=="-m": # matching probes
        onlyMethProbesWithAffyProbes = True
    elif o=="-s": # significant probes
        allRows = False

reader = csv.reader(open(infile), delimiter="\t")

writer = csv.writer(open(outfile, "w"), delimiter="\t")

###

TSS_TTS_Distance = 5000
TTS_TTS_Distance_Human = str(TSS_TTS_Distance / 1000) + "kb"

Small_TSS_TTS_Distance = 1000
Small_TTS_TTS_Distance_Human = str(Small_TSS_TTS_Distance / 1000) + "kb"

# load data

genedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")

genes = Ensembl.ReverseGeneMapping(genedata)

exons = Ensembl.ReverseExonMapping(genedata)

transcriptionSites = Ensembl.TranscriptionSites(genedata)

cpgIslands = ExtendedBed(os.path.expanduser("~/mount/publicdata/hg18/cpgislands/cpgislands-0-index.bed"))

affyannotation = NetAffxAnnotation()

paddedGenes = Ensembl.ReverseGeneMapping(genedata, tssPadding = TSS_TTS_Distance)


def isUpstream(distance,  strand):
    if strand == "+":
        return 'Y' if distance >= 0 else 'N'
    elif strand == "-":
        return'Y' if distance <= 0 else 'N'
    else:
        # wtf went wrong here
        exit(-1)
        
def isDownstream(distance,  strand):
    if strand == "+":
        return 'Y' if distance <= 0 else 'N'
    elif strand == "-":
        return'Y' if distance >= 0 else 'N'
    else:
        # wtf went wrong here
        exit(-1)


if not affyComparisonFile == None:
    #affyMapping = ExtendedBed(os.path.expanduser("~/mount/publicdata/positions2affy/HG-U133Plus2.csv"), chrPos=0, startPos = 2, stopPos=3,  defaultkeys=["chr", "strand", "start", "stop", "affy"])
    #print affyMapping.getValuesOfOverlappingIntervals("chr16", 72982016, 72983513)
    affyComparison = IndexedCSV(affyComparisonFile)

headerRow = ['Index', 'ColumnID', 'Symbol',  'Chr', 'Mapinfo','Coord']
headerRow.extend(['PD30.Avg',   'PD56.Avg', 'Fold change', 'Log2MethFC', 'Bonferroni(p-value (PD56 vs. PD30))', 'Meth'])
headerRow.extend(['In Gene',   'Genes', 'Names', 'Gene Bounds'])
headerRow.extend([TTS_TTS_Distance_Human+' up or Gene Body',   TTS_TTS_Distance_Human+' up or Gene Body Genes', TTS_TTS_Distance_Human+' up or Gene Body Names', TTS_TTS_Distance_Human+' up or Gene Body Gene Bounds'])
headerRow.extend(['Gene TSS Distance',
                                #'Gene TTS Distance',
                                'Nearest TSS',
                                'Nearest TSS Strand',
                                #'Nearest TTS',
                                #'Nearest TTS Strand'
                                ]) #'Nearest TSS Distance', 'Nearest TTS Distance'])
headerRow.extend([TTS_TTS_Distance_Human+" TSS",
                            Small_TTS_TTS_Distance_Human+" TSS",
                            #TTS_TTS_Distance_Human+" TTS",
                            #TTS_TTS_Distance_Human+" up chr of nearest TSS",
                            #TTS_TTS_Distance_Human+" down chr of nearest TSS",
                         #TTS_TTS_Distance_Human+" up chr of nearest TTS",
                         #TTS_TTS_Distance_Human+" down chr of nearest TTS"
                            ])
headerRow.extend([TTS_TTS_Distance_Human+" upstream of nearest TSS",
                               TTS_TTS_Distance_Human+" downstream of nearest TSS",
                               Small_TTS_TTS_Distance_Human+" upstream of nearest TSS",
                               Small_TTS_TTS_Distance_Human+" downstream of nearest TSS",
                            #TTS_TTS_Distance_Human+" upstream of nearest TTS",
                            #TTS_TTS_Distance_Human+" downstream of nearest TTS"
                            ])
#headerRow.extend([TTS_TTS_Distance_Human+" of TSS of Gene it's on", 
                            #TTS_TTS_Distance_Human+" of TTS of Gene it's on"
#                            ])
headerRow.extend([ 'In Exon',  'Exons', 'In Intron', 'Intergenic'])
headerRow.extend(['In CPG Island','cpg.start', 'cpg.end', 'cpg.name', 'cpg.length', 'cpg.cpgNum', 'cpg.gcNum', 'cpg.perCpg', 'cpg.perGc', 'pgp.obsExp'])

if not affyComparisonFile == None:
    headerRow.append("Meth Probe Location on Gene with Affy probe")
    headerRow.append("Ensid")
    for key in affyComparison.keys:
        headerRow.append(key)

writer.writerow(headerRow)

for row in reader:
    try:
        index = int(row[0])
        columnid = row[1]
        symbol = row[2]
        
        chr = row[3]
        mapinfo = int(row[4]) # 1 indexed
        coord = mapinfo-1 # 1 indexed convert to 0 indexed
        pd30 = float(row[5])
        pd56 = float(row[6])
        fold = float(row[7])
        pvalue = float(row[8])
        
        methProbeSignificant = True if abs(fold) >= 1.1 and pvalue <= 0.05 else False
        
        if not chr.startswith("chr"):
            chr = "chr" + chr
        
    except ValueError:
        continue # this will be the header

    # starting information
    outputRow = [str(index), columnid, symbol, chr, str(mapinfo),str(coord)]
    
    outputRow.extend([pd30, pd56, fold, math.log(abs(fold), 2) * (1 if fold > 0 else -1), pvalue])
        
    if pvalue>0.05 or abs(fold)<1.1:
        outputRow.append("NoChange")
    else:
        outputRow.append("Hypo" if fold <0.0 else "Hyper")
    
    # genes
    ingenes = genes.getValuesOfOverlappingIntervals(chr, coord, coord)
    
    if len(ingenes)==0:
        outputRow.append("N")
    else:
        outputRow.append('Y')
    
    outputRow.append(", ".join(ingenes))

    # gene names
    geneNames = []
    for gene in ingenes:
        if genedata[gene].name != "--":
            geneNames.append(genedata[gene].name)
    
    outputRow.append(", ".join(geneNames))
        
    # gene bounds
    genebounds = []
    for geneid in ingenes:
        genebounds.append(str(genedata[geneid].start)+"-"+str(genedata[geneid].end))
    outputRow.append(", ".join(genebounds))
    
    # near genes
    inPaddedGenes = paddedGenes.getValuesOfOverlappingIntervals(chr, coord, coord)
    
    if len(inPaddedGenes)==0:
        outputRow.append("N")
    else:
        outputRow.append('Y')
    
    outputRow.append(", ".join(inPaddedGenes))

    # gene names
    geneNames = []
    for gene in inPaddedGenes:
        if genedata[gene].name != "--":
            geneNames.append(genedata[gene].name)
    
    outputRow.append(", ".join(geneNames))
        
    # gene bounds
    genebounds = []
    for geneid in inPaddedGenes:
        genebounds.append(str(genedata[geneid].start)+"-"+str(genedata[geneid].end))
    outputRow.append(", ".join(genebounds))
    
    
    # TSS distances
    tsses = []
    for gene in ingenes:
        tsses.append(str(genedata[gene].tss()-coord))

    outputRow.append(", ".join(tsses))
    
    
    # TTS distances
    ttses = []
    for gene in ingenes:
        ttses.append(str(genedata[gene].tts()-coord))
#
#    outputRow.append(", ".join(ttses))

    # nearest TSS
    nearestTSS = transcriptionSites.getNearestStartSite(chr, coord)
    nearestTSSAbs = abs(nearestTSS[0].key - coord)
    tssValues = []
    tssStrand = []
    for tss in nearestTSS:
        tssValues.append(str(tss.key-coord))
        tssStrand.append(genedata[tss.data].strand)
    outputRow.append(", ".join(tssValues))
    outputRow.append(", ".join(tssStrand))

    # nearest TTS
    nearestTTS = transcriptionSites.getNearestTerminationSite(chr, coord)
    nearestTTSAbs = abs(nearestTTS[0].key - coord)
    ttsValues = []
    ttsStrand = []
    for tts in nearestTTS:
        ttsValues.append(str(tts.key-coord))
        ttsStrand.append(genedata[tss.data].strand)
#    outputRow.append(", ".join(ttsValues))
#    outputRow.append(", ".join(ttsStrand))
    
    # nearest TSS & TTS distance
    #outputRow.append(abs(nearestTSS))
    #outputRow.append(abs(nearestTTS))
    
    
    
    # TSS
    if nearestTSSAbs<=TSS_TTS_Distance:
        withindistofTSS = "Y"
        withindistUpChrmofTSS = "Y" if min(tssValues)<=0 else "N"
        withindistDownChrmofTSS = "Y" if max(tssValues)>=0 else "N"
    else:
        withindistofTSS = "N"
        withindistUpChrmofTSS = ""
        withindistDownChrmofTSS = ""
    
    # small TSS
    if nearestTSSAbs<=Small_TSS_TTS_Distance:
        smallWithindistofTSS = "Y"
        smallWithindistUpChrmofTSS = "Y" if min(tssValues)<=0 else "N"
        smallWithindistDownChrmofTSS = "Y" if max(tssValues)>=0 else "N"
    else:
        smallWithindistofTSS = "N"
        smallWithindistUpChrmofTSS = ""
        smallWithindistDownChrmofTSS = ""

    #TTS
    if nearestTTSAbs<=TSS_TTS_Distance:
        withindistofTTS = "Y"
        withindistUpChrmofTTS = "Y" if min(ttsValues)<=0 else "N"
        withindistDownChrmofTTS = "Y" if max(ttsValues)>=0 else "N"
    else:
        withindistofTTS = "N"
        withindistUpChrmofTTS = ""
        withindistDownChrmofTTS = ""

    outputRow.extend([withindistofTSS, smallWithindistofTSS,
                      #withindistofTTS,
                      #withindistUpChrmofTSS,withindistDownChrmofTSS,
                      #withindistUpChrmofTTS, withindistDownChrmofTTS
                      ])



    # up / down stream of nearest TSS
    tssUpstream = []
    tssDownstream = []
    for i in range(len(tssValues)):
        distance = int(tssValues[i])
        strand = tssStrand[i]
        if abs(distance) <= TSS_TTS_Distance:
            tssUpstream.append(isUpstream(distance, strand))
            tssDownstream.append(isDownstream(distance, strand))
        else:
            tssUpstream.append('')
            tssDownstream.append('')
    
    outputRow.extend([", ".join(tssUpstream), ", ".join(tssDownstream) ])
    
    # small up / down stream of nearest TSS
    smallTssUpstream = []
    smallTssDownstream = []
    for i in range(len(tssValues)):
        distance = int(tssValues[i])
        strand = tssStrand[i]
        if abs(distance) <= Small_TSS_TTS_Distance:
            smallTssUpstream.append(isUpstream(distance, strand))
            smallTssDownstream.append(isDownstream(distance, strand))
        else:
            smallTssUpstream.append('')
            smallTssDownstream.append('')
    
    outputRow.extend([", ".join(smallTssUpstream), ", ".join(smallTssDownstream) ])
    
#    # up / down stream of nearest TTS
#    ttsUpstream = []
#    ttsDownstream = []
#    for i in range(len(ttsValues)):
#        distance = int(ttsValues[i])
#        strand = ttsStrand[i]
#        if abs(distance) <= TSS_TTS_Distance:
#            ttsUpstream.append(isUpstream(distance, strand))
#            ttsDownstream.append(isDownstream(distance, strand))
#        else:
#            ttsUpstream.append('')
#            ttsDownstream.append('')
#    
#    outputRow.extend([", ".join(ttsUpstream), ", ".join(ttsDownstream) ])
    

    # TTS / TTS of gene it's on
    
#    nearTSSofGene = False
#    for tss in tsses:
#        if abs(int(tss))<TSS_TTS_Distance:
#            nearTSSofGene = True
#    
#    nearTTSofGene = False
#    for tts in ttses:
#        if abs(int(tts))<TSS_TTS_Distance:
#            nearTTSofGene = True
#    
#    if len(ingenes)==0:
#        # not on a gene
#        outputRow.extend(["",
#          #""
#          ])
#    else:
#        outputRow.extend(["Y" if nearTSSofGene else "N", 
#        #"Y" if nearTTSofGene else "N"
#        ])

    # exons
    inexons = exons.getValuesOfOverlappingIntervals(chr, coord, coord)
    
    if len(inexons)==0:
        outputRow.append("N")
    else:
        outputRow.append('Y')
    
    outputRow.append(", ".join(inexons))
    
    # introns
    outputRow.append("Y" if (len(ingenes)>0 and len(inexons)==0) else "N")
    
    # intergenic
    outputRow.append("Y" if (len(ingenes)==0 and nearestTSSAbs>TSS_TTS_Distance and nearestTTSAbs>TSS_TTS_Distance) else "N")
    
    
    # cpg islands
    incpg = cpgIslands.getValuesOfOverlappingIntervals(chr, coord, coord)
    
    if len(incpg)==0:
        outputRow.append("N")
    else:
        outputRow.append('Y')
    
    # be paranoid and assume it could be in multiple cpg islands (this shouldnt ever be the case but who knows what the bed file could contain
    cpg_Starts = []
    cpg_Ends = []
    cpg_Names = []
    cpg_Lengths = []
    cpg_cpgNum = []
    cpg_gcNum = []
    cpg_perCpg = []
    cpg_perGc = []
    cpg_obsExp = []
    
    for cpg in incpg:
        cpg_Starts.append(cpg['chromStart'])
        cpg_Ends.append(cpg['chromEnd'])
        cpg_Names.append(cpg['name'])
        cpg_Lengths.append(cpg['length'])
        cpg_cpgNum.append(cpg['cpgNum'])
        cpg_gcNum.append(cpg['gcNum'])
        cpg_perCpg.append(cpg['perCpg'])
        cpg_perGc.append(cpg['perGc'])
        cpg_obsExp.append(cpg['obsExp'])
    
    outputRow.append(", ".join(cpg_Starts))
    outputRow.append(", ".join(cpg_Ends))
    outputRow.append(", ".join(cpg_Names))
    outputRow.append(", ".join(cpg_Lengths))
    outputRow.append(", ".join(cpg_cpgNum))
    outputRow.append(", ".join(cpg_gcNum))
    outputRow.append(", ".join(cpg_perCpg))
    outputRow.append(", ".join(cpg_perGc))
    outputRow.append(", ".join(cpg_obsExp))
        
    if affyComparisonFile == None:
        writer.writerow(outputRow)
    else:
            inaffys = {}
            
            for inPaddedGene in inPaddedGenes:
                
                for affyprobe in affyannotation.getAffysForEnsembl(inPaddedGene):
                    inaffys[affyprobe] = inPaddedGene

            if len(inaffys) == 0:
                # no affys found that match up but we output the line if we arent only interested in affy probe match ups
                if onlyMethProbesWithAffyProbes==False:
                    writer.writerow(outputRow)
            else:
                # print out every affy
                for inaffy in inaffys: # one row per affy

                    affyrow = outputRow[:] # clone row
                    
                    # where is the meth probe in relation to the gene that this affy probe measures
                    
                    genechr = genedata[inaffys[inaffy]].chr
                    genestrand = genedata[inaffys[inaffy]].strand
                    genestart = genedata[inaffys[inaffy]].start
                    geneend = genedata[inaffys[inaffy]].end
                                        
                    assert genestrand == "+" or genestrand == "-"
                    tss = genestart if genestrand == "+" else geneend
                    if abs(tss - coord) < TSS_TTS_Distance:
                        
                        if abs(tss-coord) < Small_TSS_TTS_Distance: # within 1kb of the TSS in either direction
                            category = "Promoter"
                        elif isUpstream(tss-coord, genestrand): # is upstream (this catches all of the 5kb upstream)
                            category = "Promoter"
                        else:
                            category = "GeneBody"
                        
                    else:
                        category = "GeneBody"
                        
                    if len(incpg)>0:
                        category = category + "-cpgIsland"
                        
                    affyrow.append(category)
                    
                    affyrow.append(inaffys[inaffy])
                        
                    affycomparisonrow = affyComparison[inaffy]
                    for key in affyComparison.keys:
                        affyrow.append(affycomparisonrow[key])
                    
                    affyProbeSignificant = True if float(affycomparisonrow['BY-fdr'])<=0.05 and abs(float(affycomparisonrow['fc']))>=1.5 else False
                    
                    # only significant probes when allRows == False
                    if (methProbeSignificant == True and affyProbeSignificant == True) or allRows == True:
                        writer.writerow(affyrow)
