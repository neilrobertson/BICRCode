import csv
import os
from datastructures.genomeintervaltree import GenomeIntervalTree
from Bio import SeqIO
import subprocess
import math

# making gene and transcript dicts should allow us to do GeneDict["GeneID"]["TranscriptID"]["ExonID"] etc

class EnsemblExon():
    def __init__(self, exonid, chr, start, end):
        if not chr.startswith("chr"):
            chr = "chr" + chr
        
        self.id = exonid
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
    
    def __repr__(self):
        return self.chr+":"+str(self.start)+"-"+str(self.end)

class EnsemblTranscript(dict):
    def __init__(self, transcriptid, geneid, chr, start, end):
        if not chr.startswith("chr"):
            chr = "chr" + chr
        
        self.id = transcriptid
        self.geneid = geneid
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
    
    def addExon(self, exon):
        self[exon.id]=exon

class EnsemblGene(dict):
    def __init__(self, genedata, geneid, chr, strand, start, end, name):
        
        self.genedata = genedata
        
        if not chr.startswith("chr"):
            chr = "chr" + chr
        
        if strand == '1':
            strand = '+'
        elif strand == '-1':
            strand = "-"

        assert strand == '+' or strand == '-', strand
        
        self.id = geneid
        self.chr = chr
        self.strand = strand
        self.start = int(start)
        self.end = int(end)
        self.name = name
    
    def addTranscript(self, transcript):
        self[transcript.id]=transcript
    
    def getPromotorRegion(self, upstreamPadding = 5000,  downstreamPadding = 1000):
        assert self.strand == '+' or self.strand == '-'
        if self.strand == '+':
            start = self.start - upstreamPadding
            end = self.start + downstreamPadding
        else:
            start = self.end - downstreamPadding
            end = self.end + upstreamPadding
        
        assert end-start == (upstreamPadding + downstreamPadding)
        return max(start,0), max(end,0)
    
    def getGeneWithPromotor(self, upstreamPadding = 5000):
        (promotorstart, promotorend) = self.getPromotorRegion()
        (genestart, geneend) = self.getBodyRegion(padding=0)
        
        return (min(promotorstart, genestart), max(promotorend, geneend))
    
    def getBodyRegion(self, padding = 1000):
        assert self.strand == '+' or self.strand == '-'
        assert abs(self.end - self.start) >= padding
        if self.strand == '+':
            start = self.start + padding
            end = self.end
        else:
            start = self.start
            end = self.end - padding
        
        assert end-start == abs(self.end - self.start) - padding
        return (start, end)
    
    def tss(self):
        assert self.strand == '+' or self.strand == '-'
        if self.strand == '+':
            return self.start
        else:
            return self.end
        
    def tts(self):
        assert self.strand == '+' or self.strand == '-'
        if self.strand == '+':
            return self.end
        else:
            return self.start
    
    def consensusExons(self):
        allexons = {}
        for transcript in self:
            for exon in self[transcript]:
                if self[transcript][exon] not in allexons:
                    allexons[self[transcript][exon]] = 0
                allexons[self[transcript][exon]]+=1
        consensusExons = []
        for exon in allexons:
            if allexons[exon] == len(self):
                consensusExons.append(exon)
        
        assert self.strand == "+" or self.strand == "-",  "Unknown strand"
        
        if self.strand == "+":
            consensusExons.sort(key=lambda exon : exon.start)
        else:
            consensusExons.sort(key=lambda exon : -1*exon.end)
        
        return consensusExons
    
    def allExons(self):
        allexons = set()
        for transcript in self:
            for exon in self[transcript]:
                if self[transcript][exon] not in allexons:
                    allexons.add(exon)
        
        assert self.strand == "+" or self.strand == "-",  "Unknown strand"
        
        allexonslist = list(allexons)
        if self.strand == "+":
            allexonslist.sort(key=lambda exon : exon.start)
        else:
            allexonslist.sort(key=lambda exon : -1*exon.end)
        
        return allexonslist
    

# positionString formatted as 'X:1-100'
def splitPosition(positionString):
    chr, startstop = positionString.split(":")
    start, end = startstop.split("-")
    
    if not chr.startswith("chr"):
        chr = "chr" + chr
    start = int(start)
    end = int(end)
    
    return (chr, start, end)

# Typical input lines:
# ENSG000XX 1:1-100 -1  ENST0000YY  1:1-100 ENSE0000ZY  1:1-10
# ENSG000XX 1:1-100 -1  ENST0000YY  1:1-100 ENSE0000ZZ  1:20-30
# ENSG000XY X:1-100 1   ENST0000YZ  X:1-100 ENSE0000YY  X:1-10
# ENSG000XY X:1-100 1   ENST0000YZ  X:1-100 ENSE0000YZ  X:20-30


# Download from Biomart as TSV:

#Ensembl Gene ID
#Chromosome Name
#Gene Start (bp)
#Gene End (bp)
#Strand
#Ensembl Transcript ID
#Transcript Start (bp)
#Transcript End (bp)
#Ensembl Exon ID
#Exon Chr Start (bp)
#Exon Chr End (bp)


defaultannotations = {"hg18" : "ncbi36.1", # Ensembl 54,
                      "hg19" : "EnsemblGenes73",
                      "hg19" : "EnsemblGenes73", # Ensembl xx
                      "mm9" : "ncbi37" # Ensembl xx
                      }


class EnsemblGenes(dict):
    def __init__(self, assembly=None, annotation=None):
        
        self.assembly = assembly
        
        if annotation == None:
            annotation = defaultannotations[assembly]
        
        self.annotation = annotation
        
        #Gets to the path for the genome build annotations
        self.base = os.path.expanduser("~/mount/publicdata/"+assembly+"/"+annotation)
        
        self.transcripts = {}
        self.exons = {}
        
        #Gets gene names
        friendlyNamesLocation = self.base+"/genenames.csv"
        friendlyNames = FriendlyGeneNames(friendlyNamesLocation)
        #Gets gene IDs from gene names
        self.friendlyNamesToIDs = FriendlyNamesToGeneIDs(friendlyNamesLocation)
        
        #Gets gene structure annotations
        self.structureFileLocation = self.base+"/genes-and-exons.csv"
        
        #Gets the columns
        gfilesettings = {"geneid" : 0, "genecoord" : 1, "genestrand" : 2,       #gene
                                  "transcriptid" : 3, "transcriptcoord":4,               #transcript
                                  "exonid":5, "exoncoord":6,                                      #exon        
                                  "delimiter" : '\t' }
        
        #Opens the gene structures
        with open(self.structureFileLocation, "r") as inputFile:
            # value of 0 = input is zero indexed, value of 1 = input is one indexed
            # always converts to zero indexed
            index = 1
            rows = csv.reader(inputFile, delimiter=gfilesettings["delimiter"])
            
            for row in rows:
                if row[0].startswith("#"):
                    continue # skip header / comments
                
                geneid = row[gfilesettings["geneid"]]
                
                if geneid not in self:
                    (chr, start, end) = splitPosition(row[gfilesettings["genecoord"]])
                    if geneid in friendlyNames:
                        name = friendlyNames[geneid]
                    else:
                        name = "--"
                    self[geneid] = EnsemblGene(self, geneid, chr, row[gfilesettings["genestrand"]], start-index, end, name)
                
                gene = self[geneid]
                
                transcriptid = row[gfilesettings["transcriptid"]]
                
                if transcriptid not in gene:
                    if transcriptid not in self.transcripts:
                        (chr, start, end) = splitPosition(row[gfilesettings["transcriptcoord"]])
                        
                        self.transcripts[transcriptid] = EnsemblTranscript(transcriptid, geneid, chr, start-index, end)
                    gene.addTranscript(self.transcripts[transcriptid])
                
                transcript = gene[transcriptid]
                
                exonid = row[gfilesettings["exonid"]]
                
                # note that we keep a global list of all exons here so that if an exon is reused
                # across multiple transcripts then we use the reference to the same exon instance
                if exonid not in transcript:
                    if exonid not in self.exons:
                        (chr, start, end) = splitPosition(row[gfilesettings["exoncoord"]])
                        self.exons[exonid] = EnsemblExon(exonid, chr, start-index, end)
                    transcript.addExon(self.exons[exonid])

    def getGeneIDs(self,name):
        if name in self.friendlyNamesToIDs:
            return self.friendlyNamesToIDs[name]
        else:
            return []

# ----------

class FriendlyGeneNames(dict):
    def __init__(self, friendlyGenesFileLocation):
        self.genes = self
        with open(friendlyGenesFileLocation, "r") as inputFile:        
            genesToUseInputFile = csv.reader(inputFile, delimiter='\t')
            for row in genesToUseInputFile:
                self[row[0]]=row[1]
                
class FriendlyNamesToGeneIDs(dict):
    def __init__(self, friendlyGenesFileLocation):
        with open(friendlyGenesFileLocation, "r") as inputFile:
            genesToUseInputFile = csv.reader(inputFile, delimiter='\t')
            for row in genesToUseInputFile:
                if row[1] in self:
                    self[row[1]].append(row[0])
                else:
                    self[row[1]]=[row[0]]

# -----------
    
class ReversePromotorMapping(GenomeIntervalTree):
    def __init__(self, ensemblGenes, upstreamPadding = 0, downstreamPadding = 0):
        super(ReversePromotorMapping, self).__init__()
        
        for geneid in ensemblGenes:
                chr = ensemblGenes[geneid].chr
                (start, end) = ensemblGenes[geneid].getPromotorRegion(upstreamPadding = upstreamPadding,  downstreamPadding = downstreamPadding)
                
                self.insertInterval(chr, start, end, geneid)

# ----------

class ReverseGeneMapping(GenomeIntervalTree):
    def __init__(self, ensemblGenes, tssPadding = 0, ttsPadding = 0):
        super(ReverseGeneMapping, self).__init__()
        
        for geneid in ensemblGenes:
                
                chr = ensemblGenes[geneid].chr
                start = ensemblGenes[geneid].start
                end = ensemblGenes[geneid].end
                strand = ensemblGenes[geneid].strand
                
                assert strand == '+' or strand == '-'
                if strand == '+':
                    start -= tssPadding
                    end += ttsPadding
                else:
                    start -= ttsPadding
                    end += tssPadding
                
                self.insertInterval(chr, start, end, geneid)

# -------------

class ReverseExonMapping(GenomeIntervalTree):
    
    def __init__(self, ensemblGenes):
        super(ReverseExonMapping, self).__init__()

        for exonid in ensemblGenes.exons:
            
                chr = ensemblGenes.exons[exonid].chr
                start = ensemblGenes.exons[exonid].start
                end = ensemblGenes.exons[exonid].end
            
                self.insertInterval(chr, int(start), int(end), exonid)


# ------------


class Node(object):
    def __init__(self,key,gene):
        self.key = key
        self.data = gene

class TranscriptionSites(object):
    
    def __init__(self, geneMapping, onlyGenes = None):
        self.startSites = GenomeIntervalTree()
        self.terminationSites = GenomeIntervalTree()
        
        self.getGeneLocations(geneMapping,onlyGenes)
    
    # read in the gene id mapping
    def getGeneLocations(self, ensemblgenes, onlyGenes):
        
        if onlyGenes == None:
            genelist = ensemblgenes
        else:
            genelist = onlyGenes

        for gene in genelist:
            self.startSites.insertInterval(ensemblgenes[gene].chr, ensemblgenes[gene].tss(), ensemblgenes[gene].tss()+1, Node(ensemblgenes[gene].tss(),gene))
            self.terminationSites.insertInterval(ensemblgenes[gene].chr, ensemblgenes[gene].tts(),ensemblgenes[gene].tts()+1, Node(ensemblgenes[gene].tts(),gene))
            
    def getNearestSites(self, tree, chrm, start,stop,max):
        
        if not chrm.startswith("chr"):
            chrm = "chr" + chrm
        
        nodes = tree.getValuesInRange(chrm,start,stop)
        
        nodes.extend(tree.getValues(tree.before(chrm,start,max_dist=max)))
        nodes.extend(tree.getValues(tree.after(chrm,stop,max_dist=max)))
        
        if not len(nodes) > 0:
            return []
        
        seenKeys = set()
        uniqueNodes = []
        
        for node in nodes:
            if node.key not in seenKeys:
                seenKeys.add(node.key)
                uniqueNodes.append(node)
        
        return uniqueNodes
    
#    def getNearestAnySite(self,chrm,key):
#        nearestStart = self.getNearestStartSite(chrm, key)
#        nearestTermination = self.getNearestTerminationSite(chrm, key)
#        
#        assert len(nearestStart)>=1
#        assert len(nearestTermination)>=1
#        
#        nearestStart = nearestStart[0]
#        nearestTermination = nearestTermination[0]
#        
#        startDist = abs(nearestStart.key-key)
#        termDist = abs(nearestTermination.key-key)
#        
#        return nearestStart if startDist <= termDist else nearestTermination
    
    def getNearestStartSites(self, chrm, start,stop,max=500000000):
        return self.getNearestSites(self.startSites, chrm, start,stop,max)
        
    def getNearestTerminationSites(self, chrm, start,stop,max=500000000):
        return self.getNearestSites(self.terminationSites, chrm, start,stop.max)
    

class EnsemblProteinSequence(dict):
    def __init__(self, assembly="hg18", annotation="ncbi36.1"):
        
        self.genedata = EnsemblGenes(assembly,annotation)
        
        sequenceLocation = os.path.expanduser("~/mount/publicdata/"+assembly+"/"+annotation+"/protein_sequences.txt")

        sequences = SeqIO.parse(open(sequenceLocation,"r"),"fasta")
        
        self.transcriptmapping = dict()
        
        for sequence in sequences:
            
            name = sequence.name
            seq = sequence.seq.tostring()
            
            if len(name.split("|"))!=2:
                continue
            
            transcriptid,proteinid = name.split("|")
            
            self[proteinid] = seq
            
            self.transcriptmapping[transcriptid] = proteinid
    
    def hasSequence(self,id):
        return id in self.transcriptmapping or id in self
    
    def getProteinSequence(self,id):
        if id in self.transcriptmapping:
            return self[self.transcriptmapping[id]]
        else:
            return self[id]
    
    def getIEP(self,id):
        if not self.hasSequence(id):
            return None
        
        sequence = self.getProteinSequence(id)
        
        sequence_handle = open("sequence", "w")
        sequence_handle.write(sequence)
        sequence_handle.close()
        
        iep = subprocess.Popen(args=["iep", "-sequence", "sequence", "-outfile", "stdout"],stdout=subprocess.PIPE)
        output, errors = iep.communicate()
            
        iep = str(output).splitlines()[1].split(" = ")[1] # text munging, hacky
        
        return float(iep)
    
    def getGeneIEP(self,id):
        ieps = []
        for transcript in self.genedata[id]:
            iep = self.getIEP(transcript)
            if iep != None:
                ieps.append(iep)
        if len(ieps) == 0:
            return None
        else:
            return math.fsum(ieps)/float(len(ieps))
        
if __name__ == "__main__":
    
    genes = EnsemblGenes()
    
    genecoords = csv.writer(open(genes.base+"/all-genecoords.bed","w"),delimiter="\t",lineterminator="\n")
    
    for gene in genes:
        genecoords.writerow([genes[gene].chr, genes[gene].start, genes[gene].end])
    
    
    promotorcoords = csv.writer(open(genes.base+"/all-promotorcoords.bed","w"),delimiter="\t",lineterminator="\n")
    
    for gene in genes:
        promotorstart,promotorend = genes[gene].getPromotorRegion()
        promotorcoords.writerow([genes[gene].chr, promotorstart, promotorend])
    
    exoncoords = csv.writer(open(genes.base+"/all-exoncoords.bed","w"),delimiter="\t",lineterminator="\n")
    
    for exon in genes.exons:
        exoncoords.writerow([genes.exons[exon].chr, genes.exons[exon].start, genes.exons[exon].end])
        
    # mergeBed -i all-genecoords.bed | grep -v "_" > genic-regions.bed
    # mergeBed -i all-exoncoords.bed | grep -v "_" > exonic-regions.bed
    # mergeBed -i all-promotorcoords.bed | grep -v "_" > promotor-regions.bed