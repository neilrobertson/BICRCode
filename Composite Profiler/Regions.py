'''
Created on 25 Jul 2012

@author: mcbryan
'''

import genemapping.geneslicer as geneslicer
from genemapping.chrmEnds import ChromosomeEnds

from csvfile.genelist import GeneList
from genemapping.pointlist import Pointlist as PointList
from bed.treatment import BedFile as BedList

#####

class RegionBehaviour(object):
    
    def __init__(self):
        self.regionIterator = None
    
    # used to determine if a supplied object is valid to include
    def valid(self,obj):
        assert False, "Region validity not defined"
    
    # used to convert a list of values into an average
    def getAvg(self,values):
        assert False, "Get Average approach not defined"
    
    # used to get the chromosome of a gene object (region specific)
    def getChr(self,obj):
        assert False, "Get Chromosome not defined"
        
    # used to slice a gene object into slices (region specific)
    def slicer(self,obj):
        assert False, "Slicer not defined"  


####

from genemapping.Affy import GenesMapping as AffyGenesMapping

class AffyRegions(RegionBehaviour):
    def __init__(self,genemappingName,chunks,ud,ui,dd,di):

        self.geneNumbChunks = 20 if (chunks==None) else chunks
        self.upstreamDistance = 5000 if (ud==None) else ud
        self.upstreamInterval = 1000 if (ui==None) else ui
        self.downstreamDistance = 5000 if (dd==None) else dd
        self.downstreamInterval = 1000 if (di==None) else di
        
        self.genesmapping = AffyGenesMapping(genemappingName) # eg : "hg18toU133P2"
        
        self.regionIterator = GeneList
    
    def valid(self,gene):
        return gene in self.genesmapping
            
    def getChr(self,gene):
        (chrm, strand, start, end) = self.genesmapping[gene]
        return chrm
    
    def slicer(self,gene):
        return geneslicer.sliceGene(self.genesmapping[gene],self.geneNumbChunks, self.upstreamDistance, self.upstreamInterval, self.downstreamDistance, self.downstreamInterval)

####

from genemapping import Ensembl



class EnsemblRegions(RegionBehaviour):
    def __init__(self,build,annotation,exons,minExons,numbExons,tsses,chunks,ud,ui,dd,di):

        assert not (exons and tsses)
        
        #Stores all parameters
        self.exons = exons
        self.minExons = minExons
        self.numbExons = numbExons
        
        self.tsses = tsses
        
        self.geneNumbChunks = 20 if (chunks==None) else chunks
        self.upstreamDistance = 5000 if (ud==None) else ud
        self.upstreamInterval = 1000 if (ui==None) else ui
        self.downstreamDistance = 5000 if (dd==None) else dd
        self.downstreamInterval = 1000 if (di==None) else di
        
        # Gets dictionaries of genes, exons transcripts etc from Ensembl class
        #assert build == "hg18", "Non hg-18 not supported for Ensembl regions at the moment"
        
        #Gets the lengths of chromosomes
        # we keep this as a non instance variable as well so that the subclass below can access it
        genedata = Ensembl.EnsemblGenes(assembly=build, annotation=annotation)
        self.genedata = genedata
        
        self.chromosomeEnds = ChromosomeEnds(build)
        
        # Gets an object of the genelist class
        class EnsemblNameAndIDFormatter(list):
            def __init__(self, genesToUseLocation):
                geneList = GeneList(genesToUseLocation)
                
                self.seengenes = set()
                
                for gene in geneList:
                    if gene in self.seengenes:
                        # seen this gene id before (generally shouldnt be the case as using GeneList which guarantees this
                        # for the source list at least
                        pass
                    elif gene in genedata:
                        # we've not seen the gene but it's in the ensembl ids list
                        self.append(gene)
                        self.seengenes.add(gene)
                    else:
                        # it's not in the ensembl ids list, it could be a gene name
                        found = False
                        for geneid in genedata.getGeneIDs(gene):
                            if geneid not in self.seengenes:
                                self.append(geneid)
                                self.seengenes.add(geneid)
                                found = True
                        if not found:
                            print "No GeneID for:"+gene
                                                        
                print genesToUseLocation + ":" +str(len(self))
        
        self.regionIterator = EnsemblNameAndIDFormatter
    
    def valid(self, gene):
        #Tests if the gene in the input file is in the list of genes/features created in Regions.
        if gene not in self.genedata:
            print gene + " not in genedata"
            return False
        #Tests if the gene is the at a position which exists within the chromosome
        if self.getChr(gene) not in self.chromosomeEnds:
            return False
        return True
        
            
    def getChr(self,gene):
        return self.genedata[gene].chr
    
    def addPrefix(self,keys,slices,prefix):
        newkeys = []
        for key in keys:
            newkeys.append(prefix+key)
        for slice in slices:
            slice.sliceid = prefix+slice.sliceid
        return newkeys,slices
        
    def sliceTSS(self,gene,distance,interval):
        #if genedata[gene].strand == "+":
        return geneslicer.slicePoint(self.genedata[gene].tss(), distance, interval, distance, interval,self.genedata[gene].strand)
        #else:
        #    return geneslicer.slicePoint(genedata[gene].tss(), distance*-1, interval*-1, distance*-1, interval*-1)
    
    def sliceGene(self,gene, geneNumbChunks, upstreamDistance, upstreamInterval, downstreamDistance, downstreamInterval):
        return geneslicer.sliceGene((self.genedata[gene].chr, self.genedata[gene].strand, self.genedata[gene].start, self.genedata[gene].end),
                                    geneNumbChunks, upstreamDistance, upstreamInterval, downstreamDistance, downstreamInterval)
    
    def sliceExons(self,gene):
        # get numbExon exon slices
        exonslices = []
        sliceid = 0
        previousPoint = None
        consensusExons = self.genedata[gene].consensusExons()
        if len(consensusExons)<self.minExons:
            return [], []
        for exon in consensusExons:
            sliceid+=1
            assert self.genedata[gene].strand == "+" or self.genedata[gene].strand == "-",  "Unknown strand"
            if previousPoint != None:
                thisPoint = exon.start if self.genedata[gene].strand == "+" else exon.end
                exonslices.append(geneslicer.GeneSlice("Intron "+str(sliceid-1), min(thisPoint,previousPoint),  max(thisPoint, previousPoint)))
            previousPoint = exon.end if self.genedata[gene].strand == "+" else exon.start
            exonslices.append(geneslicer.GeneSlice("Exon "+str(sliceid), exon.start, exon.end))
        keys, slices = [slice.sliceid for slice in exonslices[:self.numbExons*2-1]], exonslices[:self.numbExons*2-1]
        return keys, slices
    
    #Takes a feature and slices it
    #Uses the up and downstream distances and intervals parameters
    #returns a gene slicer object which contains an array of keys (one for each slice), and a start and end coordinate for the feature
    def slicer(self,gene):
        # Tests if a TSS
        if self.tsses:
            keys, slices = self.sliceTSS(gene,max(self.downstreamDistance,self.upstreamDistance),min(self.downstreamInterval,self.upstreamInterval))
            keys, slices = self.addPrefix(keys,slices,"TSS-")
        # Otherwise its a gene
        else:
            #Gets lists of keys and slices, each is a string - containing later information to be split?
            keys, slices = self.sliceGene(gene,self.geneNumbChunks,self.upstreamDistance,self.upstreamInterval,self.downstreamDistance,self.downstreamInterval)
            keys, slices = self.addPrefix(keys,slices,"Gene-")

        return keys, slices
    

#####

class PointRegions(RegionBehaviour):
    def __init__(self,build,ud,ui,dd,di):

        self.upstreamDistance = 100000 if (ud==None) else ud
        self.upstreamInterval = 1000 if (ui==None) else ui
        self.downstreamDistance = 100000 if (dd==None) else dd
        self.downstreamInterval = 1000 if (di==None) else di
        
        self.chromosomeEnds = ChromosomeEnds(build)
        
        self.regionIterator = PointList
    
    def valid(self,point):
        return self.getChr(point) in self.chromosomeEnds
            
    def getChr(self,point):
        (chrm, point) = point
        return chrm
    
    def slicer(self,point):
        # print point if slices will be out of bounds
        chrm, position = point
        if (position - self.upstreamDistance < 0) or (position+self.downstreamDistance > self.chromosomeEnds[chrm]):
            print chrm + ":" + str(position - self.upstreamDistance) + "-" + str(position + self.downstreamDistance)
        
        return geneslicer.slicePoint(position, self.upstreamDistance, self.upstreamInterval, self.downstreamDistance, self.downstreamInterval)
    
    
#####

class BedRegions(RegionBehaviour):
    def __init__(self,build,chunks,ud,ui,dd,di):

        self.numbChunks = 20 if (chunks==None) else chunks
        self.upstreamDistance = 100000 if (ud==None) else ud
        self.upstreamInterval = 1000 if (ui==None) else ui
        self.downstreamDistance = 100000 if (dd==None) else dd
        self.downstreamInterval = 1000 if (di==None) else di
        
        print self.numbChunks, self.upstreamDistance, self.upstreamInterval, self.downstreamDistance, self.downstreamInterval
        
        self.chromosomeEnds = ChromosomeEnds(build)
        
        self.regionIterator = BedList
    
    def valid(self,bedentry):
        return self.getChr(bedentry) in self.chromosomeEnds
            
    def getChr(self,bedentry):
        return bedentry.chr
    
    def slicer(self,bedentry):
        return geneslicer.sliceGene([bedentry.chr,bedentry.strand,bedentry.start,bedentry.stop],self.numbChunks, self.upstreamDistance, self.upstreamInterval, self.downstreamDistance, self.downstreamInterval)


