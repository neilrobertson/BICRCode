#!/usr/bin/env python

import gc
import os
import sys
import csv
import getopt
import re
from classification import *
from bed.treatment import Bed as BedIntervalTree
from genemapping import Ensembl
from bed.treatment import SimpleBed, ExtendedBed
from csvfile.genelist import GeneList
from sequence.genome import MemoryGenome as Genome
from csvfile.indexedcsv import IndexedCSV
from csvfile.indexedcsv import ColumnIndex
from affy.NetAffxAnnotation import NetAffxAnnotation
from genemapping.chrmEnds import ChromosomeEnds
import math

gc.disable()


try:
    opts, args = getopt.getopt(sys.argv[1:], "i:o:g:r:c:e:a:", [])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    sys.exit(2)

# takes in a csv file of intervals and tells us some stuff about them
inputFileIsGenes = False

genelists = []
genelists_bound = {}
genelists_bound_promotors = {}

regions = []

controlAffyExpressionData = None
rnaSeqExpressionData = None

inGeneFilter = None

assembly = "hg18"

for o, a in opts:
    if o=="-i":
        infile = a
    elif o=="-o":
        outfile = a
    elif o=="-g":
        genelists.append(GeneList(a))
    elif o=="-r":
        if (a.count(",")>0):
            a,padding = a.split(",")
            regions.append(BedIntervalTree(a,padding=int(padding)))
        else:
            regions.append(BedIntervalTree(a))
    elif o=="-c":
        controlAffyExpressionData = IndexedCSV(a)
    #Annotated difference file input
    elif o=="-e":
        rnaSeqExpressionData = IndexedCSV(a,key="test_id")
    elif o=="-a":
        assembly = a
    
    
UPSTREAM_PROMOTOR_DIST = 2000
DOWNSTREAM_PROMOTOR_DIST = 2000

writer = csv.writer(open(outfile, "w"), delimiter="\t")

genome = Genome(assembly)

###

# load data

genedata = Ensembl.EnsemblGenes(assembly=assembly)

genes = Ensembl.ReverseGeneMapping(genedata)

genespluspromotor = Ensembl.ReverseGeneMapping(genedata, tssPadding = UPSTREAM_PROMOTOR_DIST)

genepromotors = Ensembl.ReversePromotorMapping(genedata, upstreamPadding = UPSTREAM_PROMOTOR_DIST, downstreamPadding = DOWNSTREAM_PROMOTOR_DIST)

exons = Ensembl.ReverseExonMapping(genedata)

transcriptionSites = Ensembl.TranscriptionSites(genedata)

# UCSC table browser - Expression & Regulation - CpG Islands
# Download all columns with exception of "bin"
cpgIslands = ExtendedBed(os.path.expanduser("/mnt/50tb/publicdata/"+assembly+"/CpGIslands/cpgislands.bed"))

# UCSC table browser - Mapping and Sequencing - Chromosome Bands
# Download all columns with exception of "gieStain"
gBanding = ExtendedBed(os.path.expanduser("/mnt/50tb/publicdata/"+assembly+"/G-Banding/cytogenetic.map.bed"),defaultkeys=["chrm","start","stop","band"],forcekeys=True)

chromosomeEnds = ChromosomeEnds(assembly)

###
###
###

headerRow = ['id','chr', 'start', 'stop']


intervals = SimpleBed(infile)

# extra columns from the input file
headerRow.extend(intervals.header[3:])

headerRow.extend(['In genes', "In genes + promotor",  "intergenic", 
                "iswholelyinintron", "isinpromotor","isindownstreampromotor", "isinupstreampromotor", 'overlapsTSS',
                'In exons', 'in.cpg'])
# "Ebox Motif", "Ebox Motif - Canonical A", "Ebox Motif - Canonical G", "AP1 Motif","AP-2 Motif"])

for genelist in genelists:
    headerRow.append(genelist.getFullName())
    headerRow.append(genelist.getFullName()+"-promotoronly")
    genelists_bound[genelist.getFullName()] = set()
    genelists_bound_promotors[genelist.getFullName()] = set()

for regionIntervalTree in regions:
    headerRow.append(regionIntervalTree.getFullName())
    headerRow.append("% " + regionIntervalTree.getFullName())
    

headerRow.extend(['Genes',  'Gene Names', "Genes + promotor", "Genes + promotor Names",  'Exons',  'Nearest TSS Distance', 'Region is x of TSS', 'Nearest Gene ID', 'Nearest Gene Loc'])

headerRow.extend(['cpg.start', 'cpg.end', 'cpg.name', 'cpg.length', 'cpg.cpgNum', 'cpg.gcNum', 'cpg.perCpg', 'cpg.perGc', 'cpg.obsExp'])

headerRow.extend(["GC Percent"])

headerRow.extend(["band","Nearest.Chr.End", "Distance.In"])
#headerRow.extend(["Padded G-Count",  "Padded C-Count",  "Padded A-Count", "Padded T-Count",  "Padded GC Percent",  "Padded Sequence Ebox" ])



if controlAffyExpressionData != None:
    # make a list of all the ensembl genes we have on the array
    affyannotation = NetAffxAnnotation(genome = "hg18", array="HG-U133_Plus_2", version="29")
    
    ensembls = set()
    
    for probe in controlAffyExpressionData:
        ensembls.update(affyannotation.getValues(probe, "Ensembl"))
    
    onArrayTranscriptionSites = Ensembl.TranscriptionSites(genedata,onlyGenes = ensembls)
        
    ensemblExpressionValues = dict()
    
    for ensembl in ensembls:
        affys = affyannotation.getAffysForEnsembl(ensembl)

        total = math.fsum([float(controlAffyExpressionData[affy]["Value"]) for affy in affys])
        ensemblExpressionValues[ensembl] = total / float(len(affys))
        
    headerRow.extend(["Control Expression"])

if rnaSeqExpressionData != None:
    headerRow.extend(["value_1", "value_2","Nearest Gene Diff", "Nearest Gene Fold", "Nearest Gene Sig"])


#headerRow.extend(["Sequence" ])

writer.writerow(headerRow)




id = 1

print "Starting..."

for row in intervals:
    
    chr, start, stop = row
    
    classifiee = Classifiee(chr, start, stop)
    
    extradata = intervals.getExtraData(row)
    
    if isinstance(extradata,list):
        # dont know
        pass
    elif isinstance(extradata, dict):
        for key in extradata:
            classifiee.addColumn(Column(key,extradata[key]))
    else:
        assert False, "No idea what extra data is meant to be"
    
    classifiee.addColumn(Column("id", id))
    
    id += 1
    
    size = stop - start
    
    classifiee.addColumn(Column("size", size))
    
    # genes
    ingenes = genes.getValuesInRange(chr, start, stop)
    
    if len(ingenes)==0:
        classifiee.addColumn(Column("In genes", False))
    else:
        classifiee.addColumn(Column("In genes", True))
    
    classifiee.addColumn(ListColumn("Genes",  ingenes))


    # gene names
    geneNames = []
    for gene in ingenes:
        if gene in genedata:
            geneNames.append(genedata[gene].name)
    
    classifiee.addColumn(ListColumn("Gene Names", geneNames))


    inexons = exons.getValuesInRange(chr, start, stop)
    
    if len(inexons)==0:
        classifiee.addColumn(Column("In exons", False))
    else:
        classifiee.addColumn(Column("In exons", True))
    
    classifiee.addColumn(ListColumn("Exons",  inexons))


    # nearest TSS(es)
    nearestTSSes = transcriptionSites.getNearestStartSites(chr, start, stop)
    
    # find out what the stream and distance is from interval is
    TSSdistances = []
    for tss in nearestTSSes:
        distance = min(abs(start-tss.key), abs(stop-tss.key))
        if tss.key >= start and tss.key <= stop:
            distance = 0
        
        genechr = genedata[tss.data].chr
        genestrand = genedata[tss.data].strand
        genestart = genedata[tss.data].start
        geneend = genedata[tss.data].end
        geneid = genedata[tss.data].id
        geneloc = str(genechr) + ":" + str(genestart) + "-" + str(geneend)
        
        
        # also need our position relative to it
        assert genestrand == '+' or genestrand == '-',  "Unknown strand"
        
        if start <= tss.key and stop >= tss.key:
            stream = "overlaps"
        elif start < tss.key: #dont need separate branch to check stop as by previous check it must be on the same side
            assert stop < tss.key
            stream = "upstream" if genestrand == '+' else "downstream"
        elif start > tss.key: # ditto
            assert stop > tss.key
            stream = "upstream" if genestrand == '-' else "downstream"
        else:
            assert False, "Unknown stream direction"
            
        TSSdistances.append((distance, genestrand, stream, geneid, geneloc))
    
    # find the closest one(s)
    minTSSDistance = sys.maxint
    candidateTSSes = []
    
    for distance, strand, stream, geneid, geneloc in TSSdistances:
        if distance < minTSSDistance: # its closer than anything we've found so far so we discard previous
            candidateTSSes = []
            candidateTSSes.append((distance, strand, stream, geneid, geneloc))
            minTSSDistance = distance
        elif distance == minTSSDistance: # its the same distance as the other closest ones weve found so add to previous
            candidateTSSes.append((distance, strand, stream, geneid, geneloc))
        else: # larger than the previous candidates we seen so we discard it
            continue
    
    if minTSSDistance == sys.maxint:
        print chr, str(candidateTSSes),  str(nearestTSSes)
    
    tssDistances = []
    tssStreams = set()
    tssStreamsList = []
    tssGeneIDs = []
    tssGeneLoc = []
    
    # these are all the closest TSSes we could find and should all have the same distance, but lets check anyway
    for distance, strand, stream, geneid, geneLoc in candidateTSSes:
        assert distance == minTSSDistance, "TSS Distance Not Equal to Min TSS distance"
        
        tssDistances.append(str(distance))
        tssStreamsList.append(stream)
        tssStreams.add(stream)
        tssGeneIDs.append(geneid)
        tssGeneLoc.append(geneLoc)
    
    classifiee.addColumn(Column("Nearest TSS Distance", minTSSDistance))
    classifiee.addColumn(ListColumn("Region is x of TSS", tssStreams))
    classifiee.addColumn(ListColumn("Nearest Gene ID", tssGeneIDs))
    classifiee.addColumn(Column("Nearest Gene Loc", geneLoc))
    
    #Gets the nearest gene expression data where required
    if rnaSeqExpressionData != None:
        geneExpression = ColumnIndex(rnaSeqExpressionData, "ensemblid") 
        
        if len(geneExpression[tssGeneIDs[0]]) > 0:
            # Just takes the first gene in the list where multiple exist
            for i in geneExpression[tssGeneIDs[0]]:
                classifiee.addColumn(Column("value_1",rnaSeqExpressionData[i]["value_1"]))
                classifiee.addColumn(Column("value_2",rnaSeqExpressionData[i]["value_2"]))
                classifiee.addColumn(Column("Nearest Gene Diff",rnaSeqExpressionData[i]["log2(fold_change)"]))
                classifiee.addColumn(Column("Nearest Gene Fold",float(rnaSeqExpressionData[i]["value_2"]) - float(rnaSeqExpressionData[i]["value_1"])))
                classifiee.addColumn(Column("Nearest Gene Sig",rnaSeqExpressionData[i]["significant"]))
        else:
            classifiee.addColumn(Column("value_1","NA"))
            classifiee.addColumn(Column("value_2","NA"))
            classifiee.addColumn(Column("Nearest Gene Diff","NA"))
            classifiee.addColumn(Column("Nearest Gene Fold","NA"))
            classifiee.addColumn(Column("Nearest Gene Sig","NA"))

    # genes
    ingenespluspromotors = genespluspromotor.getValuesInRange(chr, start, stop)
    
    if len(ingenespluspromotors)==0:
        classifiee.addColumn(Column("In genes + promotor", False))
    else:
        classifiee.addColumn(Column("In genes + promotor", True))
    
    classifiee.addColumn(ListColumn("Genes + promotor",  ingenespluspromotors))


    genesPlusPromotorsNames = []
    for gene in ingenespluspromotors:
        if gene in genedata:
            genesPlusPromotorsNames.append(genedata[gene].name)
    
    classifiee.addColumn(ListColumn("Genes + promotor Names", genesPlusPromotorsNames))
    


    isingene = (len(ingenes)!=0)
    isinexon = (len(inexons)!=0)
    isingeneorpromotor = (len(ingenespluspromotors)!=0)
    iswholelyinintron = isingene and not isinexon

    overlapsTSS = False

    isinpromotor = False
    for gene in ingenespluspromotors:
        if gene in genedata:
            (genepromotorstart, genepromotorend) = genedata[gene].getPromotorRegion(upstreamPadding = UPSTREAM_PROMOTOR_DIST,  downstreamPadding = DOWNSTREAM_PROMOTOR_DIST)
            # does our region overlap with this bit of the promotor
            if genepromotorstart < stop and start < genepromotorend:
                isinpromotor = True
            if genedata[gene].tss() > start and genedata[gene].tss() < stop:
                overlapsTSS = True
    
    isinupstreampromotor = False
    for gene in ingenespluspromotors:
        if gene in genedata:
            (geneupstreampromotorstart, geneupstreampromotorend) = genedata[gene].getPromotorRegion(upstreamPadding = UPSTREAM_PROMOTOR_DIST,  downstreamPadding = 0)
            # does our region overlap with this bit of the promotor
            if geneupstreampromotorstart < stop and start < geneupstreampromotorend:
                isinupstreampromotor = True
    
    if isinupstreampromotor:
        assert isinpromotor
    
    isindownstreampromotor = False
    for gene in ingenespluspromotors:
        if gene in genedata:
            (genedownstreampromotorstart, genedownstreampromotorend) = genedata[gene].getPromotorRegion(upstreamPadding = 0,  downstreamPadding = DOWNSTREAM_PROMOTOR_DIST)
            assert abs(genedownstreampromotorend-genedownstreampromotorstart)==DOWNSTREAM_PROMOTOR_DIST
            # does our region overlap with this bit of the promotor
            if genedownstreampromotorstart < stop and start < genedownstreampromotorend:
                isindownstreampromotor = True
    
    if isindownstreampromotor:
        assert isinpromotor
    
    if overlapsTSS:
        assert isinpromotor
        assert isindownstreampromotor,  str(start) + "-" + str(stop) + " -- " + str(genedownstreampromotorstart) + "-" + str(genedownstreampromotorend)
        assert isinupstreampromotor, str(start) + "-" + str(stop) + " -- " + str(geneupstreampromotorstart) + "-" + str(geneupstreampromotorend)
    
    multipledifferentpromotors = False
    if isindownstreampromotor and isinupstreampromotor and not overlapsTSS:
        multipledifferentpromotors = True
    
    intergenic = not isingene and not isinexon and not isinupstreampromotor

    if isingeneorpromotor:
        assert isingene or isinpromotor, chr + ":" + str(start) + "-" + str(stop)
    
    if isinupstreampromotor:
        assert isingeneorpromotor
    
    if isinupstreampromotor:
        assert minTSSDistance <= UPSTREAM_PROMOTOR_DIST, str(TSSdistances) + str(row) + ":" + str(minTSSDistance)
    
    if isindownstreampromotor:
        assert minTSSDistance <= DOWNSTREAM_PROMOTOR_DIST, str(TSSdistances) + str(row) + ":" + str(minTSSDistance)
        
    if isingene:
        assert isingeneorpromotor
    
    classifiee.addColumn(Column("multiple different promotors", multipledifferentpromotors))
    classifiee.addColumn(Column("intergenic", intergenic))
    classifiee.addColumn(Column("overlapsTSS", overlapsTSS))
    classifiee.addColumn(Column("iswholelyinintron", iswholelyinintron))
    classifiee.addColumn(Column("isinpromotor", isinpromotor))
    classifiee.addColumn(Column("isindownstreampromotor", isindownstreampromotor))
    classifiee.addColumn(Column("isinupstreampromotor", isinupstreampromotor))
    
    

    
    
    # on a gene list regulated?
    for genelist in genelists:
        boundneargene = False        
        # for each gene at this region
        for genepluspromotor in ingenespluspromotors:
            # for each gene in the gene list            
            for pattern in genelist:
                if re.match(pattern,genedata[genepluspromotor].name):
                    boundneargene = True
                    
        classifiee.addColumn(Column(genelist.getFullName(), boundneargene))
    
    inpromotors = genepromotors.getValuesInRange(chr, start, stop)
    
    for genelist in genelists:
        boundonpromotor = False
        for promotor in inpromotors:
            for pattern in genelist:
                if re.match(pattern,genedata[promotor].name):
                    boundonpromotor = True

        classifiee.addColumn(Column(genelist.getFullName()+"-promotoronly",boundonpromotor))
    
    
    for regionIntervalTree in regions:
        overlappingintervals = regionIntervalTree.getIntervalsInRange(chr, start, stop)

        classifiee.addColumn(Column(regionIntervalTree.getFullName(), len(overlappingintervals)>0))

        overlappingBP = 0

        for overlappinginterval in overlappingintervals:
            overlappingBP += min(overlappinginterval.end,stop) - max(overlappinginterval.start,start)
        
        percentoverlap = float(overlappingBP) / float(stop-start) * 100.0
        
        classifiee.addColumn(Column("% " + regionIntervalTree.getFullName(), percentoverlap))
        
    
    
    # cpg islands
    incpg = cpgIslands.getValuesInRange(chr, start, stop)
    
    if len(incpg)==0:
        classifiee.addColumn(Column("in.cpg", False))
    else:
        classifiee.addColumn(Column("in.cpg", True))
    
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

    
    classifiee.addColumn(ListColumn("cpg.start", cpg_Starts))
    classifiee.addColumn(ListColumn("cpg.end", cpg_Ends))
    classifiee.addColumn(ListColumn("cpg.name", cpg_Names))
    classifiee.addColumn(ListColumn("cpg.length", cpg_Lengths))
    classifiee.addColumn(ListColumn("cpg.cpgNum", cpg_cpgNum))
    classifiee.addColumn(ListColumn("cpg.gcNum", cpg_gcNum))
    classifiee.addColumn(ListColumn("cpg.perCpg", cpg_perCpg))
    classifiee.addColumn(ListColumn("cpg.perGc", cpg_perGc))
    classifiee.addColumn(ListColumn("cpg.obsExp", cpg_obsExp))

    sequence = genome.getSequence(chr, start, stop)
    
    padded_sequence = genome.getSequence(chr, max(start-100,0), stop+100)

    classifiee.addColumn(Column("Sequence", sequence))
    classifiee.addColumn(Column("100b Padded Sequence", padded_sequence))

    #G-Banding
    gBand = gBanding.getValuesInRange(chr, start, stop)
    band_name = []
    
    for band in gBand:
        band_name.append(band['band'])
        
    classifiee.addColumn(ListColumn("band", band_name))
    
    #Distance from the nearest chromosome end
    chrLength = chromosomeEnds[chr]
    if stop < chrLength/2:
        classifiee.addColumn(Column("Nearest.Chr.End", start))
    else:
        classifiee.addColumn(Column("Nearest.Chr.End", chrLength-stop))
    classifiee.addColumn(Column("Distance.In", (float(start + ((stop-start)/2))/float(chrLength))*100))
    
        

    def baseCount(seq):
        freq = {'G' : 0, 'C' : 0,  'A' : 0, 'T' : 0}
        
        for base in seq:
            ub = base.upper()
            if ub not in freq:
                freq[ub] = 0
            freq[ub] = freq[ub] + 1
        
        return {'G':freq['G'], 
                      'C':freq['C'], 
                      'A':freq['A'], 
                      'T':freq['T'],
                      }
    
    sequenceCount = baseCount(sequence)
    
#    classifiee.addColumn(Column("G-Count", sequenceCount['G']))
#    classifiee.addColumn(Column("C-Count", sequenceCount['C']))
#    classifiee.addColumn(Column("A-Count", sequenceCount['A']))
#    classifiee.addColumn(Column("T-Count", sequenceCount['T']))
    
    try:
    
        classifiee.addColumn(Column("GC Percent",  (100.0 * (sequenceCount['G'] +  sequenceCount['C'])) / 
                                                    (sequenceCount['G'] + sequenceCount['C'] + sequenceCount['A'] + sequenceCount['T'])))
    
    except ZeroDivisionError:
        print chr, start, stop
        print sequence
        sys.exit()
    
    
    # ebox regex
    
    if not Genome.sequenceHasIUPAC(sequence, 'CANNTG'):
        classifiee.addColumn(Column("Ebox Motif", False))
        classifiee.addColumn(Column("Ebox Motif - Canonical G", False))
        classifiee.addColumn(Column("Ebox Motif - Canonical A", False))
    else:
        classifiee.addColumn(Column("Ebox Motif", True))
        
        if Genome.sequenceHasIUPAC(sequence, 'CACGTG'):
            classifiee.addColumn(Column("Ebox Motif - Canonical G", True))
        else:
            classifiee.addColumn(Column("Ebox Motif - Canonical G", False))
            
        if Genome.sequenceHasIUPAC(sequence, 'CACATG'):
            classifiee.addColumn(Column("Ebox Motif - Canonical A", True))
        else:
            classifiee.addColumn(Column("Ebox Motif - Canonical A", False))
        
            
    
#    if re.search('[cC][aA]..[tT][gG]',  sequence) == None:
#        classifiee.addColumn(Column("Ebox Motif", False))
#    else:
#        if re.search('[Cc][Aa][Cc][Gg][Tt][Gg]', sequence) == None:
#            classifiee.addColumn(Column("Ebox Motif", True))
#        else:
#            classifiee.addColumn(Column("Ebox Motif", "Y - Canonical"))


    if Genome.sequenceHasIUPAC(sequence, 'TGASTCA'):
        classifiee.addColumn(Column("AP1 Motif", True))
    else:
        classifiee.addColumn(Column("AP1 Motif", False))


    if Genome.sequenceHasIUPAC(sequence, 'GCCNNNGGC'):
        classifiee.addColumn(Column("AP-2 Motif", True))
    else:
        classifiee.addColumn(Column("AP-2 Motif", False))


    if controlAffyExpressionData != None:
        # find nearest ensembl
        nearestSite = onArrayTranscriptionSites.getNearestAnySite(chr, start)
        
        classifiee.addColumn(Column("Control Expression",ensemblExpressionValues[nearestSite.data]))
        
        


    ###
    ###
    ###
    
    outputRow = []
    for column in headerRow:
        outputRow.append(classifiee.getColumnValue(column))
    
    writer.writerow(outputRow)

exit(1)
