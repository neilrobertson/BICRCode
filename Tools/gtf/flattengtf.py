'''
This is created for a specific use case involving the flattening of multiple overlapping transcripts (from the same gene) into a single transcript
created from the union of all its exons.

In this case it works for histone genes which basically have only one or two exons where the first exon is
invariably shared and the second exon is optional.

Not recommended for use on the genome as a whole without checking results are sensible.

Created on 10 Jul 2013

@author: mcbryan
'''

import getopt
import sys
from genemapping.GTFReader import GTFFile
import collections

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["gtf="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    gtf = None    
    
    
    class Exon(object):
        def __init__(self,chrm,source,feature,start,end,score,strand,frame,attributes):
            
            geneid = attributes["gene_id"]
            transcriptid = attributes["transcript_id"]
            exonnumber = transcriptid + "." + attributes["exon_number"]
                
            self.chrm = chrm
            self.start = int(start)
            self.end = int(end)
            self.strand = strand
            self.geneid = geneid
            self.transcripts = set([transcriptid])
            self.exonnumber = exonnumber
        
        def extendExon(self,exon):
            assert self.chrm == exon.chrm
            assert self.strand == exon.strand
            assert self.geneid == exon.geneid
            
            self.exonnumber = self.exonnumber + "/" + exon.exonnumber
            self.start = min(self.start,exon.start)
            self.end = max(self.end,exon.end)
            self.transcripts.update(exon.transcripts)
        
        def __str__(self):
            row = [self.chrm,
                   "dexseq_prepare_annotation.py",
                   "exonic_part",
                   str(self.start),
                   str(self.end),
                   ".",
                   self.strand,
                   "."]
            attribs = ['transcripts "' + "+".join(self.transcripts) + '"',
                       'exonic_part_number "' + self.exonnumber + '"',
                       'gene_id "'+self.geneid+'"']
            row.append("; ".join(attribs))
            return "\t".join(row)
    
    class Gene(dict):
        def __init__(self):

            self.chrm = None
            self.start = None
            self.end = None
            self.strand = None
            self.geneid = None
            
            self.seenexons = []
        
        def __setitem__(self, key, value):
            if key not in self.seenexons:
                self.seenexons.append(key)
            if self.chrm == None:
                self.chrm = value.chrm
                self.start = value.start
                self.end = value.end
                self.strand = value.strand
                self.geneid = value.geneid
            else:
                self.start = min(self.start,value.start)
                self.end = max(self.end,value.end)
            return dict.__setitem__(self, key, value)
        
        def __iter__(self):
            return iter(self.seenexons)
        
        def __str__(self):
            header = [self.chrm,
                      "dexseq_prepare_annotation.py",
                      "aggregate_gene",
                      str(self.start),
                      str(self.end),
                      ".",
                      self.strand,
                      ".",
                      'gene_id "'+self.geneid+'"']
            entry = ["\t".join(header)]
            
            
            
            for exon in self.mergeExons(self.values()):
                entry.append(str(exon))
            return "\n".join(entry)
        
        
        # merge the exons that are overlapping
        def mergeExons(self,exons):
            #Sorts the list
            exons.sort(key=lambda x: (x.chrm, x.start, x.end))
           
            mergedIntervals = []
             
            #merges the list
            previous = None
            for current in exons:
               
                if previous == None:
                    previous = current
                elif (current.end >= previous.start and current.start <= previous.end and previous.chrm == current.chrm):
                    previous.extendExon(current)
                else:
                    mergedIntervals.append(previous)
                    previous = current
                
            # done forget the last one
            if previous != None:
                mergedIntervals.append(previous)
           
            return mergedIntervals
     
     
        
    for o, a in opts:
        if o=="--gtf":
            gtf = a
    

    
    genes = collections.defaultdict(Gene)
    
    for line in GTFFile(gtf):
        
        exon = Exon(*line)
        
        if exon.exonnumber not in genes[exon.geneid]:
            genes[exon.geneid][exon.exonnumber] = exon
    
    for gene in genes:
        print genes[gene]