'''
Created on 21 Jul 2010

@author: mcbryan
'''

import sys
from sam.SamFormat import SAMFile
import getopt
#from genemapping import Ensembl
from genemapping import UCSC
import collections
import math
import csv


def getOverlappingTranscripts(samEntry):
    transcripts = collections.defaultdict(list)

    for (start,end) in samEntry.getSegments():
        # each segment
        overlappingExons = exons.getValuesOfOverlappingIntervals(samEntry.chrm, start, end)
        # find the transcripts + exons for this read
        for exon in overlappingExons:
            (transcript,exonid) = exon
            transcripts[transcript].append((exonid,(start,end)))
    return transcripts

def filterAllSegmentsInExons(samEntry,transcripts):
    # are each of the segments accounted for by exons in this transcript
    filteredTranscripts = {}
    for transcript in transcripts:

        segmentsInThisTranscript = set()
        
        for (exonid,(start,end)) in transcripts[transcript]:
            segmentsInThisTranscript.add((start,end))
        
        assert len(segmentsInThisTranscript) <= len(samEntry.getSegments())
        # transcript has 1 exon per segment
        if len(segmentsInThisTranscript)==len(samEntry.getSegments()):
            filteredTranscripts[transcript] = transcripts[transcript]
    
    return filteredTranscripts


def filterConsecutiveExons(samEntry,transcripts):
    filteredTranscripts = {}
    for transcript in transcripts:
        prevId = None
        ok = True
        for (exonid,(start,end)) in transcripts[transcript]:
            assert prevId == None or exonid-prevId >= 1, "Exon id's should be incremental"
            if prevId == None or exonid-prevId==1:
                continue
            else:
                ok = False
                break
        if ok:
            filteredTranscripts[transcript] = transcripts[transcript]
    return filteredTranscripts

def filterCorrectStrand(samEntry,transcripts):
    filteredTranscripts = {}
    for transcript in transcripts:
        if samEntry.getStrand()=="." or samEntry.getStrand() == genedata[transcript].strand:
            filteredTranscripts[transcript] = transcripts[transcript]
    return filteredTranscripts


def filterIfExactMatch(samEntry,transcripts):
    hasExactMatch = False
    for transcript in transcripts:
        if samEntry.getStrand()!="." and samEntry.getSegments() > 1:
            # do any of them have a splicing boundary exactly the same as one of the splice points
            # where others dont
            for (exonid,(start,end)) in transcripts[transcript]:
                if samEntry.getStrand()!="+" and genedata[transcript][exonid].start == start:
                    hasExactMatch = True
                    break
                elif samEntry.getStrand()!="-" and genedata[transcript][exonid].end == end:
                    hasExactMatch = True
                    break
    # we now know if it has one or more transcript that is an exact match
    # so we can filter if so (otherwise we let it slide)
    if hasExactMatch:
        matchedTranscripts = {}
        for transcript in transcripts:
            if samEntry.getStrand()!="." and samEntry.getSegments() > 1:
                # do any of them have a splicing boundary exactly the same as one of the splice points
                # where others dont
                for (exonid,(start,end)) in transcripts[transcript]:
                    if samEntry.getStrand()!="+" and genedata[transcript][exonid].start == start:
                        ## add
                        matchedTranscripts[transcript]=transcripts[transcript]
                        break
                    elif samEntry.getStrand()!="-" and genedata[transcript][exonid].end == end:
                        ## add
                        matchedTranscripts[transcript]=transcripts[transcript]
                        break
        return matchedTranscripts
    else:
        return transcripts



def filterEntirelyContainedWithinExons(samEntry,transcripts):
    # look for entirely contained
    entirelyContained = False
    for transcript in transcripts:
        if samEntry.getStrand()!="." and samEntry.getSegments() > 1:                        
            # are any of them entirely contained within the exons of a transcript
            containedSoFar = True
            for (exonid,(start,end)) in transcripts[transcript]:
                if not (start >= genedata[transcript][exonid].start and end<=genedata[transcript][exonid].end):
                    containedSoFar = False
            if containedSoFar == True:
                entirelyContained = True
                break
    # now we know if there are any entirely contained or not, so we can filter based on that
    if entirelyContained:
        containedTranscripts = {}
        for transcript in transcripts:
            if samEntry.getStrand()!="." and samEntry.getSegments() > 1:                        
                # are any of them entirely contained within the exons of a transcript
                containedSoFar = True
                for (exonid,(start,end)) in transcripts[transcript]:
                    if not (start >= genedata[transcript][exonid].start and end<=genedata[transcript][exonid].end):
                        containedSoFar = False
                if containedSoFar == True:
                    containedTranscripts[transcript] = transcripts[transcript]
        return containedTranscripts
    else:
        return transcripts


def filterByPercentageCoverage(samEntry,transcripts,cutoff):
        # code for measuring coverage of read - discard transcripts with low agreement (< cutoff)
        coveredTranscripts = set()
        for transcript in transcripts:
            numbBases = 0
            coveredBases = 0
            for (exonid,(start,end)) in transcripts[transcript]:
                numbBases+=(end-start)
                exonstart = genedata[transcript][exonid].start
                exonend = genedata[transcript][exonid].end
                assert min(end,exonend) - max(start,exonstart) > 0
                coveredBases += min(end,exonend) - max(start,exonstart)
            assert coveredBases<=numbBases#                    
            percentage = float(coveredBases)/float(numbBases)
            if percentage>cutoff:
                coveredTranscripts[transcript] = transcripts[transcript]
        return coveredTranscripts




def whereIsARead(samEntry):
    
    #
    # first we look for "acceptable" assignments
    #
    overlappingTranscripts = getOverlappingTranscripts(samEntry)
    
    # each of the segments is in an exon within the transcript
    possibleTranscripts = filterAllSegmentsInExons(samEntry,overlappingTranscripts)

    # are all the exons "consecutive"
    consecutiveTranscripts = filterConsecutiveExons(samEntry,possibleTranscripts)
    
    # look for correct strand exons
    correctStrandTranscripts = filterCorrectStrand(samEntry,consecutiveTranscripts)
    
    #
    # now we start looking for "best" assignments
    #
    # exact match on an exon junction
    matchedTranscripts = filterIfExactMatch(samEntry,correctStrandTranscripts)
    
    # fully contained within the exons
    coveredTranscripts = filterAllSegmentsInExons(samEntry, matchedTranscripts)
    
    return coveredTranscripts
    
    


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", [])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
    #genedata = Ensembl.EnsemblGenes(assembly="hg18", annotation="ncbi36.1")    
    #exons = Ensembl.ReverseExonMapping(genedata)
    
    genedata = UCSC.UCSCTranscripts(assembly="hg18")
    exons = UCSC.ReverseExonMapping(genedata)
    
    
    for arg in args:
        
        # this could potentially be quite a big list for some transcripts and we might want to optimise it later
        transcriptReads = collections.defaultdict(collections.deque) # a list of the read values (fractional reads) for each transcript
        clusterReads = collections.defaultdict(collections.deque) # a list of read values for a gene (fractional reads)
        
        # count of how many segments reads have
        segmentCounts = collections.defaultdict(int)
        
        readsTotal = 0
        
        transcriptsPerRead = collections.defaultdict(int) # defaults to 0
        
        samFile = SAMFile(arg)
        # do the unique hits first
        for samEntry in samFile:
            # unique reads only first
            if not samEntry.isHeader():
                # this is a unique mapping so we can assign to an transcript without any doubt
                # so lets see where it goes
                readsTotal += 1
                segmentCounts[len(samEntry.getSegments())]+=1
                
                if readsTotal % 1000 == 0:
                    print ".",
                    if readsTotal % 10000 == 0:
                        print readsTotal
                
                transcripts = whereIsARead(samEntry)
                
                if len(transcripts) == 0:
                    # read from unknown (novel?) splicing?
                    # we probably want to save these and then try to build exons out of them later
                    pass
                
                # just some tallying for the end
                
                numbTranscripts = len(transcripts)                                
                transcriptsPerRead[numbTranscripts]+=1
                
                if numbTranscripts>0:
                    transcriptWeight = 1.0/(samEntry.getNumberHits()*numbTranscripts)
                else:
                    # read isn't anywhere
                    continue
                
                # add count to transcripts
                for transcript in transcripts:
                    if numbTranscripts == 1 and samEntry.getNumberHits() == 1:
                        # we specifically append '1.0' here so we can do float compare later
                        transcriptReads[transcript].append(1.0)
                    else:
                        transcriptReads[transcript].append(transcriptWeight)
                
                # add count to genes
                
                # get list of all cluster id's involved
                clusters = collections.defaultdict(int) # defaults to 0
                for transcript in transcripts:
                    clusters[genedata[transcript].clusterid]+=1
                
                numbClusters = len(clusters)
                
                # clusterids is now a count of how many transcripts we found have that clusterid
                for cluster in clusters:
                    if numbClusters == 1  and samEntry.getNumberHits() == 1:
                        # again this is just so that float compare works later
                        clusterReads[cluster].append(1.0)
                    else:
                        clusterReads[cluster].append(transcriptWeight*clusters[cluster])
        
        # print summaries
        
        print str(readsTotal)+" reads looked at"
        
        print transcriptsPerRead
        print segmentCounts
        
        with open(arg+"-transcripts.csv","w") as outputFile:
        
            csvout = csv.writer(outputFile,delimiter="\t")
            csvout.writerow(["TranscriptId","Symbol","AssignedCount","MinCount","MaxCount"])
        
            for transcript in transcriptReads:
                min = 0
                for read in transcriptReads[transcript]:
                    if read == 1.0:
                        min += 1 
                assigned = math.fsum(transcriptReads[transcript])
                max = len(transcriptReads[transcript])
                
                assert min <= assigned <= max
                
                csvout.writerow([transcript,genedata[transcript].geneSymbol, assigned, min, max])
                
                
        with open(arg+"-clusters.csv","w") as outputFile:

            csvout = csv.writer(outputFile,delimiter="\t")
            csvout.writerow(["ClusterId","Symbols","AssignedCount","MinCount","MaxCount"])

            for cluster in clusterReads:
                min = 0
                for read in clusterReads[cluster]:
                    if read == 1.0:
                        min += 1
                assigned = math.fsum(clusterReads[cluster])
                max = len(clusterReads[cluster])
                
                assert min <= assigned <= max
                
                csvout.writerow([cluster, ",".join(genedata.clusters[cluster].geneSymbols), assigned, min, max])
    
    print "Done"