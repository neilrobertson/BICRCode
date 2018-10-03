import os
from csvfile.indexedcsv import IndexedCSV
from datastructures.genomeintervaltree import GenomeIntervalTree
import collections

class UCSCExon():
    def __init__(self,id,chr,strand,start,end):
        self.id = id
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand
        
    def __repr__(self):
        return "("+str(self.id)+"/"+self.chr+":"+str(self.start)+"-"+str(self.end)+")"
        
class UCSCTranscript(dict):
    def __init__(self,id,chr,start,end,strand,clusterid,exonStarts,exonEnds,geneSymbol,ensemblTranscript):
        self.id = id
        self.chr = chr
        self.start = start
        self.end = end
        assert strand == "+" or strand == "-"
        self.strand = strand
        self.clusterid = clusterid
        self.geneSymbol = geneSymbol
        self.ensemblTranscript = ensemblTranscript
        
        for i in range(len(exonStarts)):
            self[i] = UCSCExon(i,chr,strand,int(exonStarts[i]),int(exonEnds[i]))
    
    def __repr__(self):
        
        exons = ""
        for exon in self:
            exons += ", " +str(self[exon])
        
        return self.id+" ("+self.geneSymbol+") @ "+self.chr+":"+str(self.start)+"-"+str(self.end) + exons

class UCSCCluster():
    def __init__(self,clusterid):
        self.clusterid = clusterid
        self.geneSymbols = set()
        self.transcripts = set()

    def addToCluster(self,transcript):
        assert transcript.clusterid == self.clusterid
        
        self.transcripts.add(transcript.id)
        self.geneSymbols.add(transcript.geneSymbol)
    
    def __str__(self):
        return ",".join(self.geneSymbols)

class UCSCTranscripts(dict):
    def __init__(self,assembly="hg18"):
        
        self.clusters = {}
        
        baseLocation = os.path.expanduser("~/mount/publicdata/"+assembly+"/ucsc/")
        
        xrefs = IndexedCSV(baseLocation + "ucsc-genes-xref.csv")
        clusters = IndexedCSV(baseLocation + "ucsc-genes-knownIsoforms.csv", keyPos=1)
        ensemblMappings = IndexedCSV(baseLocation + "ucsc-knownToEnsembl.csv")
        
        self.reverseEnsemblMappings = collections.defaultdict(list)
        
        transcripts = IndexedCSV(baseLocation + "ucsc-knowngenes.csv")
        
        for id in transcripts:
            
            chr = transcripts[id]["chrom"]
            
            if chr not in ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                           "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                           "chr21","chr22",
                           "chrX","chrY"]:
                continue
            
            strand = transcripts[id]["strand"]
            assert strand == "+" or strand == "-"
            
            if id in ensemblMappings:
                ensemblTranscript = ensemblMappings[id]["value"]
            else:
                ensemblTranscript = None
            
            start = int(transcripts[id]["txStart"])
            end = int(transcripts[id]["txEnd"])
            
            exonCount = int(transcripts[id]["exonCount"])
            exonStarts = transcripts[id]["exonStarts"].split(",") # comma seperated
            exonEnds = transcripts[id]["exonEnds"].split(",") # comma seperated
            # ucsc terminates with a comma which leaves a blank bit at the end after split
            exonStarts.remove('')
            exonEnds.remove('')
            
            # make sure it's consistent
            assert len(exonStarts) == len(exonEnds) and exonCount == len(exonStarts), str(transcripts[id]) 
            
            # not using these at the moment
            
            #if "proteinID" in transcripts[id]: # not always there
            #    proteinId = transcripts[id]["proteinID"]
            #alignId = transcripts[id]["alignID"]
            
            if id in xrefs:
                genesymbol = xrefs[id]["geneSymbol"]
            else:
                genesymbol = None
            
            assert id in clusters
            clusterid = clusters[id]["clusterId"]
            if clusterid not in self.clusters:
                self.clusters[clusterid] = UCSCCluster(clusterid)
            
            self[id]=UCSCTranscript(id,chr,start,end,strand,clusterid,exonStarts,exonEnds,genesymbol,ensemblTranscript)
            
            self.reverseEnsemblMappings[ensemblTranscript].append(id)
            
            self.clusters[clusterid].addToCluster(self[id])

    def getTranscriptsForEnsembl(self,ensemblid):
        return self.reverseEnsemblMappings[ensemblid]


class ReverseExonMapping(GenomeIntervalTree):
    
    def __init__(self, ucscTranscripts):
        super(ReverseExonMapping, self).__init__()

        for transcriptid in ucscTranscripts:
            for exonid in ucscTranscripts[transcriptid]:
                chr = ucscTranscripts[transcriptid][exonid].chr
                start = ucscTranscripts[transcriptid][exonid].start
                end = ucscTranscripts[transcriptid][exonid].end
            
                self.insertInterval(chr, int(start), int(end), (transcriptid, exonid))

if __name__ == "__main__":
    test = UCSCTranscripts()
    
    print test["uc003swh.1"]