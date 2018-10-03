#!/usr/bin/env python

import os
import sys
import csv
import getopt
import math
from genemapping import Ensembl
from bed.treatment import ExtendedBed
from csvfile.indexedcsv import IndexedCSV
from affy.NetAffxAnnotation import NetAffxAnnotation
from classification import *
from sequence.genome import Genome
from genemapping.chrmEnds import ChromosomeEnds
import collections


try:
    opts, args = getopt.getopt(sys.argv[1:], "i:o:a:",       [  "input=",
                                                                "output=",
                                                                "affyarray=",
                                                                
                                                                "bed=",
                                                                "bedgraphfc=",
                                                                "bedgraphbeta=",
                                                                "bedgraphp=",
                                                                "genebygene=",
                                                                "affyvalues="])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    sys.exit(2)

# takes in a csv file of (point) coordinates and tells us some stuff about them



affyComparisonFile = None
onlyMethProbesWithAffyProbes = False
allRows = True
bedTrackLoc = None
bedGraphTrackLoc = None
bedGraphPTrackLoc = None
bedGraphBetaTrackLoc = None
bedGraphFcTrackLoc = None
geneByGeneOutputLoc = None
affyValueColumns = None

for o, a in opts:
    if o in ["-i", "--input"]:
        infile = a
    elif o in ["-o", "--output"]:
        outfile = a
    elif o in ["-a","--affyarray"]:
        affyComparisonFile = a
    elif o in ["--bed"]:
        bedTrackLoc = a
    elif o in ["--bedgraphfc"]:
        bedGraphFcTrackLoc = a
    elif o in ["--bedgraphbeta"]:
        bedGraphBetaTrackLoc = a
    elif o in ["--bedgraphp"]:
        bedGraphPTrackLoc = a
    elif o in ["--genebygene"]:
        geneByGeneOutputLoc = a
    elif o in ["--affyvalues"]:
        affyValueColumns = a.split(",")

reader = IndexedCSV(infile)

writer = csv.writer(open(outfile, "w"), delimiter="\t")

###

def distanceHumanReadable(dist):
    return str(dist / 1000) + "kb"


TSS_TTS_Distance = 1000
SURROUNDING_SEQUENCE_Distance = 250 # each side
WINDOW_SIZE = 500
WINDOW_OFFSET = 5

# load data

genome = Genome(genomeBuild = "hg18")

chromosomeEnds = ChromosomeEnds("hg18")

genedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")

genes = Ensembl.ReverseGeneMapping(genedata)

exons = Ensembl.ReverseExonMapping(genedata)

transcriptionSites = Ensembl.TranscriptionSites(genedata)

cpgIslands = ExtendedBed(os.path.expanduser("~/mount/publicdata/hg18/cpgislands/cpgislands-0-index.bed"))

affyannotation = NetAffxAnnotation()

paddedGenes = Ensembl.ReverseGeneMapping(genedata, tssPadding = TSS_TTS_Distance)


# store full mapping here for the end
geneToMethProbeMapping = collections.defaultdict(list)
geneToAffyProbeMapping = collections.defaultdict(list)

def isUpstream(distance,  strand):
    assert strand in ['+','-']
    if strand == "+":
        return distance >= 0
    else:
        return distance <= 0
        
def isDownstream(distance,  strand):
    assert strand in ['+','-']
    if strand == "+":
        return distance <= 0
    else:
        return distance >= 0


if not affyComparisonFile == None:
    affyComparison = IndexedCSV(affyComparisonFile)
    def isAffySignificant(affyid):
        return float(affyComparison[affyid]['BH-fdr'])<=0.05 and abs(float(affyComparison[affyid]['fc']))>=1.5
    def affyDirection(affyid):
        return "Up" if float(affyComparison[affyid]['fc']) > 0.0 else "Down"


def sequenceGCandCpGcontent(sequence):
    
    assert len(sequence) > 0
    
    Ns = sequence.count('N')
    Gs = sequence.count('G')
    Cs = sequence.count('C')
    
    cpgs = sequence.count('CG')
    
    effectiveSeqLen = len(sequence) - Ns
    
    if effectiveSeqLen == 0: # all N's
        gcContent = 0.0
        cpgContent = 0.0
        cpgRatio = 0.0
    else:
        gcContent = float(Gs+Cs) / float(effectiveSeqLen)
        
        cpgContent = float(cpgs) / float(effectiveSeqLen)
        
        if cpgs>0:
            cpgRatio = (float(cpgs) * effectiveSeqLen) / float(Gs*Cs) # from Schubeler
        else:
            cpgRatio = 0.0 # avoid div by 0 error if Gs+Cs == 0, shortcut cpgs == 0 then cpgRatio = 0
        
#    if cpgRatio > 2.0:
#        print "---"
#        print sequence
#        print "Effective length:",effectiveSeqLen
#        print "Ns:",Ns
#        print "Gs:",Gs
#        print "Cs:",Cs
#        print "CpGs:",cpgs
#        print "CpG ratio:",cpgRatio

    return Gs, Cs, cpgs, effectiveSeqLen, gcContent, cpgContent, cpgRatio



headerRow = reader.keys[:]

headerRow.extend(['Methylation',"0 index Coord"])

# probe classification
headerRow.extend(['Gs','Cs','CpGs','Effective Seq Len','GC content',"CpG content","CpG ratio","Probe Classification"])
                                  
# in genes
headerRow.extend(['In Gene',
                  'In Genes',
                  'In Gene Names',
                  'In Gene Bounds',
                  'In Gene TSS Distances'])

# in padded genes
headerRow.extend(['In '+distanceHumanReadable(TSS_TTS_Distance)+' up or Gene Body',
                  'In '+distanceHumanReadable(TSS_TTS_Distance)+' up or Gene Body Genes',
                  'In '+distanceHumanReadable(TSS_TTS_Distance)+' up or Gene Body Names',
                  'In '+distanceHumanReadable(TSS_TTS_Distance)+' up or Gene Body Gene Bounds',
                  'In '+distanceHumanReadable(TSS_TTS_Distance)+' up or Gene Body TSS Distances'])



headerRow.extend([  'Nearest TSS',
                    'Nearest TSS Direction',
                    'Nearest TSS Distance'])


headerRow.extend([  "Upstream of nearest TSS",
                    "Downstream of nearest TSS"])


headerRow.extend(["Within "+distanceHumanReadable(TSS_TTS_Distance)+" of a TSS",
                  "Within "+distanceHumanReadable(TSS_TTS_Distance)+" upstream of a TSS",
                  "Within "+distanceHumanReadable(TSS_TTS_Distance)+" downstream of a TSS"])


headerRow.extend([ 'In Exon',  'Exons', 'In Intron', 'Intergenic'])


headerRow.extend(['In CPG Island',
                  'cpg.start',
                  'cpg.end',
                  'cpg.name',
                  'cpg.length',
                  'cpg.cpgNum',
                  'cpg.gcNum',
                  'cpg.perCpg',
                  'cpg.perGc',
                  'pgp.obsExp'])





# is there an affy array in the comparison?

if affyComparisonFile == None:
    writer.writerow(headerRow)
else:
    headerRowPlusAffy = headerRow[:]
    headerRowPlusAffy.append("Meth Probe Location on Gene with Affy probe")
    headerRowPlusAffy.append("Ensid")
    headerRowPlusAffy.extend(affyComparison.keys)
    writer.writerow(headerRowPlusAffy)



# track files
if bedTrackLoc != None:
    bedTrackFile = open(bedTrackLoc,"w")
    bedTrack = csv.writer(bedTrackFile,delimiter="\t")

if bedGraphFcTrackLoc != None:
    bedGraphFcTrackFile = open(bedGraphFcTrackLoc,"w")
    bedGraphFcTrack = csv.writer(bedGraphFcTrackFile,delimiter="\t")

if bedGraphBetaTrackLoc != None:
    bedGraphBetaTrackFile = open(bedGraphBetaTrackLoc,"w")
    bedGraphBetaTrack = csv.writer(bedGraphBetaTrackFile,delimiter="\t")
    
    
    
if bedGraphPTrackLoc != None:
    bedGraphPTrackFile = open(bedGraphPTrackLoc,"w")
    bedGraphPTrack = csv.writer(bedGraphPTrackFile,delimiter="\t")











for index in reader:
    
    chr = reader[index]["CHR"]
    if not chr.startswith("chr"):
        chr = "chr" + chr
    
    mapinfo = int(reader[index]["MAPINFO"]) # 1 indexed
    coord = mapinfo-1 # 1 indexed convert to 0 indexed
    

    classifiee = Classifiee(chr, coord, coord)
    
    classifiee.addColumn(Column("0 index Coord",str(coord)))
    
    classifiee.addColumn(Column("CHR",chr))
    
    for key in reader.keys:
        classifiee.addColumn(Column(key, reader[index][key]))
    
    
    fold = float(float(reader[index]["FC"]))
    deltabeta = float(float(reader[index]["Delta Beta"]))
    
    fdrpvalue = float(float(reader[index]["BY-fdr(tt)"]))
    
    
    logFC = math.log(abs(fold),2) * (1 if fold > 0 else -1)
    
    
    methProbeSignificant = abs(fold) >= 1.5 and abs(deltabeta) >= 0.1 and fdrpvalue <= 0.05
    classifiee.addColumn(Column("Significant", methProbeSignificant))
    
    
    
    if methProbeSignificant and fold > 0:
        direction = "Hyper"
    elif methProbeSignificant and fold < 0:
        direction = "Hypo"
    else:
        direction = "NoChange"
    classifiee.addColumn(Column("Methylation", direction))

    assert coord <= chromosomeEnds[chr] and coord >= 0


    # notice cast to uppercase
    sequence = genome.getSequence(chr,max(coord-SURROUNDING_SEQUENCE_Distance,0),
                                      min(coord+SURROUNDING_SEQUENCE_Distance,chromosomeEnds[chr])).upper()

    
    Gs, Cs, cpgs, effectiveSeqLen, gcContent, cpgContent, cpgRatio = sequenceGCandCpGcontent(sequence)

    classifiee.addColumn(Column("Gs", Gs))
    classifiee.addColumn(Column("Cs", Cs))
    classifiee.addColumn(Column("CpGs", cpgs))
    classifiee.addColumn(Column("Effective Seq Len", effectiveSeqLen))

    classifiee.addColumn(Column("GC content", gcContent))
    classifiee.addColumn(Column("CpG content", cpgContent))
    classifiee.addColumn(Column("CpG ratio", cpgRatio))
    
    
    hcp = False
    icp = False
    if gcContent >= 0.55 and cpgRatio >= 0.75:
        hcp = True
    elif cpgRatio >= 0.48:
        icp = True
        
    classifiee.addColumn(Column("Probe Classification","hcp" if hcp else ("icp" if icp else "lcp")))


    ####
    # BED TRACK
    ####
    if bedTrackLoc != None:
        # bed format used here:        
        #chrom
        #chromStart
        #chromEnd
        #name
        #score
        #strand - Defines the strand - either '+' or '-'.
        #thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
        #thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
        #itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.

        if direction == "Hyper":
            bedColour = "0,255,0"
        elif direction == "Hypo":
            bedColour = "255,0,0"
        else:
            bedColour = "0,0,0"
            
        bedTrack.writerow([chr,coord,coord+1,reader[index]["TargetID"],1000,'+',coord,coord+1,bedColour])


    if bedGraphFcTrackLoc != None:
        bedGraphFcTrack.writerow([chr,coord,coord+1,logFC])
    if bedGraphPTrackLoc != None:
        bedGraphPTrack.writerow([chr,coord,coord+1,-10.0*math.log(fdrpvalue,10)])
    if bedGraphBetaTrackLoc != None:
        bedGraphBetaTrack.writerow([chr,coord,coord+1,deltabeta])
         

    ####
    # GENES
    ####

    
    # genes
    ingenes = genes.getValuesOfOverlappingIntervals(chr, coord, coord)
        
    classifiee.addColumn(Column('In Gene', len(ingenes)>0))
    classifiee.addColumn(Column('In Genes', ", ".join(ingenes)))

    # gene names
    geneNames = []
    for gene in ingenes:
        if genedata[gene].name != "--":
            geneNames.append(genedata[gene].name)
    classifiee.addColumn(Column('In Gene Names', ", ".join(geneNames)))
        
    # gene bounds
    genebounds = []
    for geneid in ingenes:
        genebounds.append(str(genedata[geneid].start)+"-"+str(genedata[geneid].end) + "("+str(genedata[geneid].strand+")"))
    classifiee.addColumn(Column('In Gene Bounds', ", ".join(genebounds)))
    
    
    ####
    # PADDED GENES
    ####
    
    # in padded genes
    inPaddedGenes = paddedGenes.getValuesOfOverlappingIntervals(chr, coord, coord)
    
    classifiee.addColumn(Column('In '+distanceHumanReadable(TSS_TTS_Distance)+' up or Gene Body', len(inPaddedGenes)>0))
    classifiee.addColumn(Column('In '+distanceHumanReadable(TSS_TTS_Distance)+' up or Gene Body Genes', ", ".join(inPaddedGenes)))

    # padded gene names
    geneNames = []
    for gene in inPaddedGenes:
        if genedata[gene].name != "--":
            geneNames.append(genedata[gene].name)
    classifiee.addColumn(Column('In '+distanceHumanReadable(TSS_TTS_Distance)+' up or Gene Body Names', ", ".join(geneNames)))
    
    # padded gene bounds
    genebounds = []
    for geneid in inPaddedGenes:
        genebounds.append(str(genedata[geneid].start)+"-"+str(genedata[geneid].end))
    classifiee.addColumn(Column('In '+distanceHumanReadable(TSS_TTS_Distance)+' up or Gene Body Gene Bounds', ", ".join(genebounds)))
    

    for gene in inPaddedGenes:
        geneToMethProbeMapping[gene].append(classifiee)
        
    
    
    ####
    # TSSes
    ####
    
    
    # TSS distances
    tsses = []
    for gene in ingenes:
        tsses.append(str(genedata[gene].tss()-coord))
    classifiee.addColumn(Column('In Gene TSS Distances', ", ".join(tsses)))
    
    
    
    # Padded TSS distances
    paddedtsses = []
    for gene in inPaddedGenes:
        paddedtsses.append(str(genedata[gene].tss()-coord))
    classifiee.addColumn(Column('In '+distanceHumanReadable(TSS_TTS_Distance)+' up or Gene Body TSS Distances', ", ".join(paddedtsses)))
    

    

    # nearest TSS
    nearestTSS = transcriptionSites.getNearestStartSite(chr, coord)
    nearestTSSAbs = abs(nearestTSS[0].key - coord)
    tssValues = []
    tssStrand = []
    for tss in nearestTSS:
        tssValues.append(str(tss.key-coord))
        tssStrand.append(genedata[tss.data].strand)
    classifiee.addColumn(Column('Nearest TSS', ", ".join(tssValues)))
    classifiee.addColumn(Column('Nearest TSS Direction', ", ".join(tssStrand)))
    
    # nearest TSS distance
    classifiee.addColumn(Column('Nearest TSS Distance', abs(nearestTSSAbs)))
    
    # up / down stream of nearest TSS
    tssUpstream = False
    tssDownstream = False
    for i in range(len(tssValues)):
        distance = int(tssValues[i])
        strand = tssStrand[i]
        if isUpstream(distance, strand):
            tssUpstream = True
        if isDownstream(distance, strand):
            tssDownstream = True
    
    classifiee.addColumn(Column("Upstream of nearest TSS", tssUpstream))
    classifiee.addColumn(Column("Downstream of nearest TSS", tssDownstream))
    
    
    classifiee.addColumn(Column("Within "+distanceHumanReadable(TSS_TTS_Distance)+" of a TSS", nearestTSSAbs<=TSS_TTS_Distance))
    
    # up / down stream of nearest TSS within 1kb
    tssUpstream = False
    tssDownstream = False
    for i in range(len(tssValues)):
        distance = int(tssValues[i])
        strand = tssStrand[i]
        if abs(distance) <= TSS_TTS_Distance:
            if isUpstream(distance, strand):
                tssUpstream = True
            if isDownstream(distance, strand):
                tssDownstream = True

            
    classifiee.addColumn(Column("Within "+distanceHumanReadable(TSS_TTS_Distance)+" upstream of a TSS", tssUpstream))
    classifiee.addColumn(Column("Within "+distanceHumanReadable(TSS_TTS_Distance)+" downstream of a TSS", tssDownstream))

    # exons
    inexons = exons.getValuesOfOverlappingIntervals(chr, coord, coord)
    
    classifiee.addColumn(Column('In Exon',len(inexons)>0))
    classifiee.addColumn(Column('Exons',", ".join(inexons)))
    classifiee.addColumn(Column('In Intron',len(ingenes)>0 and len(inexons)==0))
    classifiee.addColumn(Column('Intergenic',len(ingenes)==0 and nearestTSSAbs>TSS_TTS_Distance))

    



    
    # cpg islands
    incpg = cpgIslands.getValuesOfOverlappingIntervals(chr, coord, coord)
    classifiee.addColumn(Column('In CPG Island',len(incpg)>0))
        
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
    
    classifiee.addColumn(Column('cpg.start',", ".join(cpg_Starts)))
    classifiee.addColumn(Column('cpg.end',", ".join(cpg_Ends)))
    classifiee.addColumn(Column('cpg.name',", ".join(cpg_Names)))
    classifiee.addColumn(Column('cpg.length',", ".join(cpg_Lengths)))
    classifiee.addColumn(Column('cpg.cpgNum',", ".join(cpg_cpgNum)))
    classifiee.addColumn(Column('cpg.gcNum',", ".join(cpg_gcNum)))
    classifiee.addColumn(Column('cpg.perCpg',", ".join(cpg_perCpg)))
    classifiee.addColumn(Column('cpg.perGc',", ".join(cpg_perGc)))
    classifiee.addColumn(Column('pgp.obsExp',", ".join(cpg_obsExp)))
        
    
    outputRow = []
    for column in headerRow:
        outputRow.append(classifiee.getColumnValue(column))
    
    
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
            for affy in inaffys: # one row per affy

                affyrow = outputRow[:] # clone row
                
                # where is the meth probe in relation to the gene that this affy probe measures
                
                genechr = genedata[inaffys[affy]].chr
                genestrand = genedata[inaffys[affy]].strand
                genestart = genedata[inaffys[affy]].start
                geneend = genedata[inaffys[affy]].end
                                    
                assert genestrand == "+" or genestrand == "-"
                tss = genestart if genestrand == "+" else geneend
                if abs(tss - coord) < TSS_TTS_Distance:
                    
                    if abs(tss-coord) < TSS_TTS_Distance: # within 1kb of the TSS in either direction
                        category = "Promoter"
                    else:
                        category = "GeneBody"
                    
                else:
                    category = "GeneBody"
                    
                if len(incpg)>0:
                    category = category + "-cpgIsland"
                    
                affyrow.append(category)
                
                affyrow.append(inaffys[affy])
                    
                affycomparisonrow = affyComparison[affy]
                for key in affyComparison.keys:
                    affyrow.append(affycomparisonrow[key])
                
                affyProbeSignificant = isAffySignificant(affy)
                
                # only significant probes when allRows == False
                if (methProbeSignificant == True and affyProbeSignificant == True) or allRows == True:
                    writer.writerow(affyrow)








def classifyDirection(d,upLabel, downLabel, noChangeLabel = "NoChange", threshold = 3):
    if d[noChangeLabel]+d[upLabel]+d[downLabel]==0:
        return "NoProbes"
    elif d[upLabel] == d[downLabel] == 0:
        return noChangeLabel
    elif d[upLabel]*threshold <= d[downLabel]: # implicit 0*threshold < 1+
        return downLabel
    elif d[downLabel]*threshold <= d[upLabel]:
        return upLabel
    else:
        return "Ambiguous"

if geneByGeneOutputLoc != None:

    # gene by gene basis
    print "Starting gene by gene"

    with open(geneByGeneOutputLoc,"w") as geneByGeneFile:
        geneByGene = csv.writer(geneByGeneFile,delimiter="\t")
        
        header = [ "EnsemblID",
                  "All NoChange","All Hyper","All Hypo", "All",
                  "Promotor NoChange","Promotor Hyper","Promotor Hypo", "Promotor",
                  "GeneBody NoChange","GeneBody Hyper","GeneBody Hypo", "GeneBody",
                  "Promotor Classification","Max CpG Window Ratio","Avg CpG Window Ratio",
                  "CpG island in Promotor", "CpG island in Gene body",
                  "Expr NoChange", "Expr Up", "Expr Down", 
                  "Expression"]
        
        chr = genedata[gene].chr
        
        # average affy value columns
        if affyValueColumns != None:
            header.extend(affyValueColumns)
        
        geneByGene.writerow(header)
        
        for gene in geneToMethProbeMapping:
            
            row = [gene]
            
            # meth probes
            methDirections =         {"NoChange":0,"Hyper":0,"Hypo":0}
            promotorMethDirections = {"NoChange":0,"Hyper":0,"Hypo":0}
            geneBodyMethDirections = {"NoChange":0,"Hyper":0,"Hypo":0}
            
            for probeClassification in geneToMethProbeMapping[gene]:
                probeMeth = probeClassification.getRawColumnValue("Methylation")
                methDirections[probeMeth]+=1
                
                withinPromotor = abs(probeClassification.start - genedata[gene].tss()) <= TSS_TTS_Distance
                
                if withinPromotor:
                    promotorMethDirections[probeMeth]+=1
                else:
                    geneBodyMethDirections[probeMeth]+=1
                
            row.extend([methDirections["NoChange"],
                        methDirections["Hyper"],
                        methDirections["Hypo"],
                        classifyDirection(methDirections,upLabel="Hyper",downLabel="Hypo")])
            
            row.extend([promotorMethDirections["NoChange"],
                        promotorMethDirections["Hyper"],
                        promotorMethDirections["Hypo"],
                        classifyDirection(promotorMethDirections,upLabel="Hyper",downLabel="Hypo")])
        
            row.extend([geneBodyMethDirections["NoChange"],
                        geneBodyMethDirections["Hyper"],
                        geneBodyMethDirections["Hypo"],
                        classifyDirection(geneBodyMethDirections,upLabel="Hyper",downLabel="Hypo")])
        
        
            # classify promotor
            
            hcp = False
            icp = False
            maxcpgratio = 0.0
            cpgratios = []
            
            # start from TSS position - TSS distance
            # go to TSS position + TSS distance + 1, but stop 1 window size away
            # do it in steps of window offset
            for leftpos in range(genedata[gene].tss()-TSS_TTS_Distance,
                                    genedata[gene].tss()+TSS_TTS_Distance+1-(WINDOW_SIZE),
                                    WINDOW_OFFSET):
                rightpos = leftpos + WINDOW_SIZE
                
                #print leftpos, rightpos, chromosomeEnds[genedata[gene].chr]
                
                # notice cast to uppercase
                sequence = genome.getSequence(genedata[gene].chr,
                                                  max(leftpos,0),
                                                  min(rightpos,chromosomeEnds[genedata[gene].chr])).upper()
                
                Gs, Cs, cpgs, effectiveSeqLen, gcContent, cpgContent, cpgRatio = sequenceGCandCpGcontent(sequence)
                
                if gcContent >= 0.55 and cpgRatio >= 0.75:
                    hcp = True
                elif cpgRatio >= 0.48:
                    icp = True
                
                maxcpgratio = max(maxcpgratio,cpgRatio)
                cpgratios.append(cpgRatio)
            
            promotorClassification = "hcp" if hcp else ("icp" if icp else "lcp")
            
            row.append(promotorClassification)
        
            row.append(maxcpgratio)
            
            avgcpgratio = math.fsum(cpgratios) / float(len(cpgratios))
        
            row.append(avgcpgratio)
        
        
        
        
            # cpg islands
            promotorincpg = cpgIslands.getValuesOfOverlappingIntervals(genedata[gene].chr, genedata[gene].tss()-TSS_TTS_Distance, genedata[gene].tss()+TSS_TTS_Distance)
            row.append(len(promotorincpg)>0)
        
        
            genebodyincpg = cpgIslands.getValuesOfOverlappingIntervals(genedata[gene].chr,
                                                                       min(genedata[gene].tss(), genedata[gene].tts()),
                                                                       max(genedata[gene].tss(), genedata[gene].tts()))
            row.append(len(genebodyincpg)>0)
        
        
        
        
        
            affyDirections = {"NoChange":0,"Up":0,"Down":0}
            
            # average affy values for this gene
            affyValues = collections.defaultdict(list)
            
            for affy in affyannotation.getAffysForEnsembl(gene):
                
                for affyValueColumn in affyValueColumns:
                    affyValues[affyValueColumn].append(float(affyComparison[affy][affyValueColumn]))
                
                if isAffySignificant(affy):
                    affyDirections[affyDirection(affy)]+=1
                else:
                    affyDirections["NoChange"]+=1
            
            row.extend([affyDirections["NoChange"],
                        affyDirections["Up"],
                        affyDirections["Down"],
                        classifyDirection(affyDirections,upLabel="Up",downLabel="Down")])
            
            # append the average affy value for each column
            row.extend([
                        (math.fsum(affyValues[affyValueColumn]) / float(len(affyValues[affyValueColumn])) if len(affyValues[affyValueColumn]) else "--") 
                        for affyValueColumn in affyValueColumns
                        ])
            
            
            
            # work out the list of windows around the promotor
            
            
            
            
            
            geneByGene.writerow(row)