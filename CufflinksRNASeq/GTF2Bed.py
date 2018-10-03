'''
Created on 18 Mar 2011

@author: mcbryan
'''

import sys
import getopt
import csv
import collections
from csvfile.indexedcsv import IndexedCSV
from filesystem.mkdir import makeDirectory

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "G:i:g:o:", [
                                                      # command args go here
                                                      "gtf=",
                                                      "isoformexp=",
                                                      "geneexp=",
                                                      "outputfolder="
                                                      ])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        print "Usage: xyz"
        sys.exit(2)
        
    isoformFile = None
    geneFile = None    
    
    for o, a in opts:
        if (o=="-G") or (o=="--gtf"):
            gtfFileLoc = a
        elif (o=="-i") or (o=="--isoformexp"):
            isoformFile = IndexedCSV(a,keyPos=0)
        elif (o=="-g") or (o=="--geneexp"):
            geneFile = IndexedCSV(a,keyPos=0)
        elif (o=="-o") or (o=="--outputfolder"):
            outputFilePrefix = a
        else:
            print "Unknown parameter: "+o+" "+a
            sys.exit(2)
    
    makeDirectory(outputFilePrefix)
    
    gtfReader = csv.reader(open(gtfFileLoc,"r"),delimiter="\t")
    
    class Transcript(list):
        def __init__(self):
            self.chr = None
            self.strand = None
        
        def appendExon(self,chr,strand,exon):
            if self.chr == None:
                self.chr = chr
            else:
                assert self.chr == chr
            
            assert strand == "+" or strand == "-" or strand == ".", strand
            if self.strand == None:
                self.strand = strand
            else:
                assert self.strand == strand
            
            self.append(exon)
    
        def getTranscriptBounds(self):
            # sort by exon number
            self.sort(key=lambda x:x[0])
            
            transcript_start = self[0][1]
            transcript_end = self[len(transcripts[transcriptid])-1][2]
            
            assert transcript_end >= transcript_start
            return transcript_start,transcript_end
    
    
    transcripts = collections.defaultdict(Transcript)
    genes = collections.defaultdict(list)
    transcript_to_genes = {}

    # load up all the data
    for row in gtfReader:
        detailsColumn = row[8]
        details = dict(item.replace("\"","").split(" ") for item in detailsColumn.split("; "))
        
        chr = row[0]
        start = int(row[3])-1 # convert to 0 base
        end = int(row[4])
        strand = row[6]

        geneid = details["gene_id"]
        transcriptid = details["transcript_id"]
        transcript_to_genes[transcriptid] = geneid
        
        genes[geneid].append(transcriptid)
        
        exon_number = int(details["exon_number"])

        transcripts[transcriptid].appendExon(chr,strand,(exon_number,start,end))
    
    with open(outputFilePrefix+"/transcripts.bed","w") as outputFile:
        outcsv = csv.writer(outputFile,delimiter="\t")
        for transcriptid in transcripts:
            
            transcript_start, transcript_end = transcripts[transcriptid].getTranscriptBounds()
            
            thickStart = transcript_start
            thickEnd = transcript_end
            
            # default colour and score
            rgb = "0,0,0"
            score = 500
            
            # colour appropriately
            if isoformFile != None and transcriptid in isoformFile and isoformFile[transcriptid]["significant"]=="yes":
                
                logtype = "ln(fold_change)" if "ln(fold_change)" in isoformFile[transcriptid] else "log2(fold_change)"
                
                rgb = "0,255,0" if float(isoformFile[transcriptid][logtype])>0 else "255,0,0"
                score = 1000
            
            row = [transcripts[transcriptid].chr,transcript_start,transcript_end,
                   transcriptid,score,transcripts[transcriptid].chr,thickStart,thickEnd,rgb,len(transcripts[transcriptid])]
            
            blockSizes = []
            blockStarts = []
            
            for (exon_number,start,end) in transcripts[transcriptid]:
                blockSizes.append(str(end-start))
                blockStarts.append(str(start-transcript_start))
            
            row.extend([",".join(blockSizes),",".join(blockStarts)])
            
            outcsv.writerow(row)
    
    
    with open(outputFilePrefix+"/transcripts.geneids.bed","w") as outputFile:
        outcsv = csv.writer(outputFile,delimiter="\t")
        for transcriptid in transcripts:
            
            transcript_start, transcript_end = transcripts[transcriptid].getTranscriptBounds()
            
            thickStart = transcript_start
            thickEnd = transcript_end
            
            # default colour and score
            rgb = "0,0,0"
            score = 500
            
            # colour appropriately
            
            if geneFile != None and geneFile[geneid]["significant"]=="yes":
                    
                logtype = "ln(fold_change)" if "ln(fold_change)" in geneFile[geneid] else "log2(fold_change)"

                rgb = "0,255,0" if float(geneFile[geneid][logtype])>0 else "255,0,0"
                score = 1000
            
            row = [transcripts[transcriptid].chr,transcript_start,transcript_end,
                   transcript_to_genes[transcriptid],score,transcripts[transcriptid].chr,thickStart,thickEnd,rgb,len(transcripts[transcriptid])]
            
            blockSizes = []
            blockStarts = []
            
            for (exon_number,start,end) in transcripts[transcriptid]:
                blockSizes.append(str(end-start))
                blockStarts.append(str(start-transcript_start))
            
            row.extend([",".join(blockSizes),",".join(blockStarts)])
            
            outcsv.writerow(row)
    
    
    if geneFile != None:
        with open(outputFilePrefix+"/genes.bed","w") as outputFile:
            outcsv = csv.writer(outputFile,delimiter="\t")
        
            for geneid in geneFile:
                
                chr,loc = geneFile[geneid]["locus"].split(":")
                
                #start,stop = loc.split("-")
                
                rgb = "0,0,0"
                score = 500
                
                if geneFile[geneid]["significant"]=="yes":
                    
                    logtype = "ln(fold_change)" if "ln(fold_change)" in geneFile[geneid] else "log2(fold_change)"
    
                    rgb = "0,255,0" if float(geneFile[geneid][logtype])>0 else "255,0,0"
                    score = 1000
                
                start = None
                stop = None
                
                for transcriptid in genes[geneid]:
                    transcript_start, transcript_stop = transcripts[transcriptid].getTranscriptBounds()
                    
                    if start == None or start>transcript_start:
                        start = transcript_start
                    if stop == None or stop<transcript_stop:
                        stop = transcript_stop
                
                row = [chr,start,stop,geneid,score,".",start,stop,rgb]
                
                outcsv.writerow(row)
    
    if isoformFile != None:
        with open(outputFilePrefix+"/transcripts.html","w") as outputFile:
            print >> outputFile, '<html><head><title>Transcripts</title><body>'
            for id in isoformFile:
                print >> outputFile, '<h2 id="'+id+'">'+id+'</h2>'
                for key in isoformFile[id]:
                    print >> outputFile, '<p><b>'+key+':</b> '+isoformFile[id][key]+'</p>'
            print >> outputFile, '</body></html>'

    if geneFile != None:
        with open(outputFilePrefix+"/genes.html","w") as outputFile:
            print >> outputFile, '<html><head><title>Genes</title><body>'
            for id in geneFile:
                print >> outputFile, '<h2 id="'+id+'">'+id+'</h2>'
                for key in geneFile[id]:
                    print >> outputFile, '<p><b>'+key+':</b> '+geneFile[id][key]+'</p>'
            print >> outputFile, '</body></html>'

    
            
    print "Done"