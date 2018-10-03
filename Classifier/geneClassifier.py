#!/usr/bin/env python

import os
import sys
import csv
import getopt
import math
import re
from classification import *

sys.path.append(os.path.expanduser("~/mount/repository/shared/baselib"))

from mapping import ensToAffyMapping
from genemapping import Ensembl
from bed.treatment import SimpleBed, Bed, ExtendedBed
from genemapping.transcriptionsites import TranscriptionSites
from csvfile.indexedcsv import IndexedCSV
from sequence.genome import Genome, UnknownChromosomeException

try:
    opts, args = getopt.getopt(sys.argv[1:], "ei:o:", [])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    sys.exit(2)

# takes in a csv file of intervals and tells us some stuff about them
for o, a in opts:
    if o=="-i":
        infile = a
    elif o=="-o":
        outfile = a

promotorUp = 2000
promotorDown = 2000

intervals = csv.reader(open(infile, "r"), delimiter="\t")

writer = csv.writer(open(outfile, "w"), delimiter="\t")

cpgIslands = ExtendedBed(os.path.expanduser("~/mount/publicdata/hg18/cpgislands/cpgislands-0-index.bed"))

genome = Genome()

# load gene data
genedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")

headerRow = [ 'Ensembl','Name','chr', 'start', 'stop', 'strand', 'No. Transcripts', 'Avg. Exons per Transcript', "Unique Exons per Gene", "Start positions", "Start positions / No. Transcripts"]

#"Promotor G-Count",  "Promotor C-Count",  "Promotor A-Count", "Promotor T-Count" , 
headerRow.extend(["Promotor region start", "Promotor region end",  "Promotor GC Percent",  ])
headerRow.extend(["Promotor Ebox","TATA within 50b upstream", "Initiator", "DPE within 50b downstream", "BREu within 60b upstream", "BREd within 60b upstream", "BRE both within 60b upstream",  "TATA and BREu"])
headerRow.extend(["Sequence"])

headerRow.extend([
                "has.cpg", 
                "cpg.start", 
                "cpg.end", 
                "cpg.name", 
                "cpg.length", 
                "cpg.cpgNum", 
                "cpg.gcNum", 
                "cpg.perCpg", 
                "cpg.perGc", 
                "cpg.obsExp"])

headerRow.extend([
                "promotor.cpg", 
                "promotor.cpg.start", 
                "promotor.cpg.end", 
                "promotor.cpg.name", 
                "promotor.cpg.length", 
                "promotor.cpg.cpgNum", 
                "promotor.cpg.gcNum", 
                "promotor.cpg.perCpg", 
                "promotor.cpg.perGc", 
                "promotor.cpg.obsExp"])

writer.writerow(headerRow)

for row in intervals:
    geneid = row[0]
        
    if geneid not in genedata:
        print "Missing Gene details:"+geneid
        continue
    
    chr = genedata[geneid].chr
    genestart = genedata[geneid].start
    genestop = genedata[geneid].end
    strand = genedata[geneid].strand
    
    classifiee = Classifiee(chr, genestart, genestop)
    
    classifiee.addColumn(Column("strand", strand))
    
    #(chr, strand, genestart, genestop) = genedata[geneid]
    
    if chr == "chrM":
        continue
    
    classifiee.addColumn(Column("Name", genedata[geneid].name))
    
    classifiee.addColumn(Column("Ensembl", geneid))
    
    classifiee.addColumn(Column("No. Transcripts", len(genedata[geneid])))
    
    exonsInTranscripts = 0
    uniqueExons = set()
    for transcriptid in genedata[geneid]:
        exonsInTranscripts += len(genedata[geneid][transcriptid])
        for exonid in genedata[geneid][transcriptid]:
            uniqueExons.add(exonid)
    
    startpos = set()
    for transcriptid in genedata[geneid]:
        # start of this transcript
        transcriptstart = genedata[geneid][transcriptid].start if strand == '+' else genedata[geneid][transcriptid].end
        startpos.add(transcriptstart)
    
    classifiee.addColumn(Column("Start positions", len(startpos)))
    
    classifiee.addColumn(Column("Start positions / No. Transcripts", float(len(startpos))/len(genedata[geneid])))
    
    classifiee.addColumn(Column("Avg. Exons per Transcript",float(exonsInTranscripts)/len(genedata[geneid])))
    
    classifiee.addColumn(Column("Unique Exons per Gene", len(uniqueExons) ))
    
    (promotorstart, promotorstop) = genedata[geneid].getPromotorRegion(upstreamPadding = promotorUp,  downstreamPadding = promotorDown)
    
    classifiee.addColumn(Column( "Promotor region start", promotorstart))
    classifiee.addColumn(Column("Promotor region end", promotorstop))
    
    try:
        sequence = genome.getSequence(chr, promotorstart, promotorstop)
    except UnknownChromosomeException:
        continue

    classifiee.addColumn(Column("Sequence", sequence))

    def baseCount(seq):
        freq = {}
        
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
    
    classifiee.addColumn(Column("Promotor G-Count", sequenceCount['G']))
    classifiee.addColumn(Column("Promotor C-Count", sequenceCount['C']))
    classifiee.addColumn(Column("Promotor A-Count", sequenceCount['A']))
    classifiee.addColumn(Column("Promotor T-Count", sequenceCount['T']))
    
    classifiee.addColumn(Column("Promotor GC Percent",  (100.0 * (sequenceCount['G'] +  sequenceCount['C'])) / 
                                                                                                    (sequenceCount['G'] + sequenceCount['C'] + sequenceCount['A'] + sequenceCount['T'])))
    
    
    
    
    # ebox regex
    if Genome.sequenceHasIUPAC(genome.getSequence(chr,genedata[geneid].getPromotorRegion(upstreamPadding = promotorUp,  downstreamPadding = promotorDown)), 'CANNTG'):
        if Genome.sequenceHasIUPAC(genome.getSequence(chr,genedata[geneid].getPromotorRegion(upstreamPadding = promotorUp,  downstreamPadding = promotorDown)), 'CACGTG'):
            classifiee.addColumn(Column("Promotor Ebox", "Y - Canonical"))
        else:
            classifiee.addColumn(Column("Promotor Ebox", "Y"))
    else:
        classifiee.addColumn(Column("Promotor Ebox", "N"))
    
    
    # tata box regex
    if Genome.sequenceHasIUPAC(genome.getSequence(chr,genedata[geneid].getPromotorRegion(upstreamPadding = 50,  downstreamPadding = 0)), 'TATAA',  strand=strand):
        classifiee.addColumn(Column("TATA within 50b upstream", "Y"))
    else:
        classifiee.addColumn(Column("TATA within 50b upstream", "N"))
    
    # initiator
    if Genome.sequenceHasIUPAC(genome.getSequence(chr,genedata[geneid].getPromotorRegion(upstreamPadding = 5,  downstreamPadding = 10)), 'YYANWYY',  strand=strand):
        classifiee.addColumn(Column("Initiator", "Y"))
    else:
        classifiee.addColumn(Column("Initiator", "N"))
    
    # downstream promotor
    if Genome.sequenceHasIUPAC(genome.getSequence(chr,genedata[geneid].getPromotorRegion(upstreamPadding = 0,  downstreamPadding = 50)), 'RGWYV',  strand=strand):
        classifiee.addColumn(Column("DPE within 50b downstream", "Y"))
    else:
        classifiee.addColumn(Column("DPE within 50b downstream", "N"))
    
    # BRE
    if Genome.sequenceHasIUPAC(genome.getSequence(chr,genedata[geneid].getPromotorRegion(upstreamPadding = 60,  downstreamPadding = 0)), 'SSRCGCC',  strand=strand):
        classifiee.addColumn(Column("BREu within 60b upstream", "Y"))
    else:
        classifiee.addColumn(Column("BREu within 60b upstream", "N"))
    
    # BRE 2
    if Genome.sequenceHasIUPAC(genome.getSequence(chr,genedata[geneid].getPromotorRegion(upstreamPadding = 60,  downstreamPadding = 0)), 'RTDKKKK',  strand=strand):
        classifiee.addColumn(Column("BREd within 60b upstream", "Y"))
    else:
        classifiee.addColumn(Column("BREd within 60b upstream", "N"))
    
    # both BRE
    if Genome.sequenceHasIUPAC(genome.getSequence(chr,genedata[geneid].getPromotorRegion(upstreamPadding = 60,  downstreamPadding = 0)), 'SSRCGCC',  strand=strand):
        if Genome.sequenceHasIUPAC(genome.getSequence(chr,genedata[geneid].getPromotorRegion(upstreamPadding = 60,  downstreamPadding = 0)), 'RTDKKKK',  strand=strand):
            classifiee.addColumn(Column("BRE both within 60b upstream", "Y"))
        else:
            classifiee.addColumn(Column("BRE both within 60b upstream", "N"))
    else:
        classifiee.addColumn(Column("BRE both within 60b upstream", "N"))
        
    
    # TATA and BRu
    if Genome.sequenceHasIUPAC(genome.getSequence(chr,genedata[geneid].getPromotorRegion(upstreamPadding = 50,  downstreamPadding = 0)), 'TATAA',  strand=strand):
        if Genome.sequenceHasIUPAC(genome.getSequence(chr,genedata[geneid].getPromotorRegion(upstreamPadding = 60,  downstreamPadding = 0)), 'SSRCGCC',  strand=strand):
            classifiee.addColumn(Column("TATA and BREu", "Y"))
        else:
            classifiee.addColumn(Column("TATA and BREu", "N"))
    else:
        classifiee.addColumn(Column("TATA and BREu", "N"))
    
    
    
    # cpg island within promotor
    incpg = cpgIslands.getValuesOfOverlappingIntervals(chr, promotorstart, promotorstop)
    
    if len(incpg)==0:
        classifiee.addColumn(Column("promotor.cpg", 'N'))
    else:
        classifiee.addColumn(Column("promotor.cpg", 'Y'))
    
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

    classifiee.addColumn(ListColumn("promotor.cpg.start", cpg_Starts))
    classifiee.addColumn(ListColumn("promotor.cpg.end", cpg_Ends))
    classifiee.addColumn(ListColumn("promotor.cpg.name", cpg_Names))
    classifiee.addColumn(ListColumn("promotor.cpg.length", cpg_Lengths))
    classifiee.addColumn(ListColumn("promotor.cpg.cpgNum", cpg_cpgNum))
    classifiee.addColumn(ListColumn("promotor.cpg.gcNum", cpg_gcNum))
    classifiee.addColumn(ListColumn("promotor.cpg.perCpg", cpg_perCpg))
    classifiee.addColumn(ListColumn("promotor.cpg.perGc", cpg_perGc))
    classifiee.addColumn(ListColumn("promotor.cpg.obsExp", cpg_obsExp))
    
    
    # cpg island within gene body
    (bodystart,  bodyend) = genedata[geneid].getBodyRegion(padding = 0)
    incpg = cpgIslands.getValuesOfOverlappingIntervals(chr, bodystart, bodyend)
    
    if len(incpg)==0:
        classifiee.addColumn(Column("has.cpg", 'N'))
    else:
        classifiee.addColumn(Column("has.cpg", 'Y'))
    
    # be paranoid and assume it could be in multiple cpg islands
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
    
    
    # ------
    
    outputRow = []
    for column in headerRow:
        outputRow.append(classifiee.getColumnValue(column))
    
    writer.writerow(outputRow)
