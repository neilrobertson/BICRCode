from datastructures.genomeintervaltree import GenomeIntervalTree
from GTFReader import GTFFile
from operator import attrgetter

# making gene and transcript dicts should allow us to do GeneDict["GeneID"]["TranscriptID"]["ExonID"] etc
class EnsemblExon():
    def __init__(self, transcript, chrm, start, end, strand):
        if not chrm.startswith("chr"):
            chrm = "chr" + chrm
        
        self.transcript = transcript
        self.chrm = chrm
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        assert self.strand in ['+','-']
        assert self.start <= self.end
    
    def __repr__(self):
        return self.chrm+":"+str(self.start)+"-"+str(self.end)

class EnsemblTranscript(list):
    def __init__(self, gene, transcriptid, transcriptname, source):
        
        self.gene = gene
        self.id = transcriptid
        self.name = transcriptname
        self.source = source
        self.strand = None
        self.chrm = None
    
    def addExon(self, exon):
        # verify strand is same as parent transcript
        assert self.strand == None or self.strand == exon.strand
        self.strand = exon.strand
        # very if chromosome is same as parent transcript
        assert self.chrm == None or self.chrm == exon.chrm
        self.chrm = exon.chrm
        
        # check the gene is consistent as well
        assert self.gene.strand == None or self.gene.strand == self.strand
        self.gene.strand = self.strand
        assert self.gene.chrm == None or self.gene.chrm == self.chrm
        self.gene.chrm = self.chrm
        
        self.append(exon)
        
        self.sortExons(reverse=True if self.strand=="-" else False)
    
    def sortExons(self,reverse=False):
        assert self.strand != None
        
        if not reverse:
            # ascending order by start coordinate
            self.sort(key=attrgetter('start', 'end'))
        else:
            # descending order by end coordinate
            self.sort(key=attrgetter('end', 'start'), reverse=True)
    
    def getTSS(self):
        return self[0].start if self.strand=="+" else self[0].end
    
    def getTTS(self):
        return self[-1].end if self.strand=="+" else self[-1].start

    def getPromotorRegion(self, upstreamPadding = 5000,  downstreamPadding = 1000):
        assert self.strand == '+' or self.strand == '-'
        if self.strand == '+':
            start = self.start - upstreamPadding
            end = self.start + downstreamPadding
        else:
            start = self.end - downstreamPadding
            end = self.end + upstreamPadding
        
        assert end-start == (upstreamPadding + downstreamPadding)
        return (start, end)

    def getLocation(self):
        return self.chrm+":"+str(min(self.getTSS(),self.getTTS())) + "-" + str(max(self.getTSS(),self.getTTS()))

    def getSource(self):
        return self.source
    
    def __repr__(self):
        return self.id + " " + self.getLocation()

class EnsemblGene(dict):
    def __init__(self, genedata, geneid, genename):
        
        self.genedata = genedata
        
        self.id = geneid
        self.name = genename
        
        self.chrm = None
        self.strand = None
        
        self.source = None
        
    def addTranscript(self, transcript):
        self[transcript.id]=transcript
        
        if self.source == None:
            self.source = transcript.source
        else:
            assert self.source == transcript.source
    
    def getPromotorRegion(self, upstreamPadding = 5000,  downstreamPadding = 1000):
        assert False, "Not implemented"
    
    def getGeneWithPromotor(self, upstreamPadding = 5000):
        assert False, "Not implemented"
    
    def getBodyRegion(self, padding = 1000):
        assert False, "Not implemented"
    
    def tss(self):
        assert False, "Not implemented"
        
    def tts(self):
        assert False, "Not implemented"
    
    def consensusExons(self):
        assert False, "Not implemented"
        
    def __repr__(self):
        return self.id + ": " + self.values().__repr__()
    



class EnsemblGenes(dict):
    def __init__(self, assembly="hg18", annotation="NCBI36.54"):
        
        self.transcripts = {}
        
        index = 1
        
        for record in GTFFile("~/mount/publicdata/"+assembly+"/rna-seq/ensembl/Homo_sapiens."+annotation+".fixed.gtf"):
            chrm,source,feature,start,end,score,strand,frame,attributes = record
            
            geneid = attributes['gene_id']
            genename = attributes['gene_name']
            
            if geneid not in self:
                self[geneid] = EnsemblGene(self, geneid, genename)
                
            gene = self[geneid]
            
            transcriptid = attributes['transcript_id']
            transcriptname = attributes['transcript_name']
            
            if transcriptid not in self.transcripts:
                self.transcripts[transcriptid] = EnsemblTranscript(gene, transcriptid, transcriptname, source)
            
            if transcriptid not in gene:
                gene.addTranscript(self.transcripts[transcriptid])
                
            transcript = gene[transcriptid]

            if feature == "exon":
                transcript.addExon(EnsemblExon(transcript, chrm, start-index, end, strand))




if __name__ == "__main__":
    
    genes = EnsemblGenes()
    
    print len(genes)
    print len(genes.transcripts)
    
    print genes['ENSG00000148187']
    
    #for gene in genes:
    #    print genes[gene]
