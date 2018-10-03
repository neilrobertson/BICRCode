# class for looking up information about a gene
# so far includes:
#	- coordinate lookup (what coord for gene x?)
#	- name lookup (given gene id, provide common name?)
#	- interval overlap (given interval, what genes overlap?)

from genemapping.coordlookup import *
from datastructures.genomicintervalfinder import *

class GeneLookup:
	def __init__(self):
		self._coordlookup = None
		self._overlaplookup = None
		self._poverlaplookup = None
		
	def initCoordsFromFile(self, filename, filesettings):
		self._coordlookup = CoordLookup()
		self._coordlookup.setFromFile(filename, filesettings)
		
	def initOverlapsFromCoords(self):
		print "converting gene coordinate lookup into overlap lookup"
		assert not self._coordlookup is None
		self._overlaplookup = GenomicIntervalFinder("Gene")
		self._overlaplookup.initFromCoordLookup(self._coordlookup)

	def initPromoterOverlapsFromCoords(self, pbounds):
		print "converting gene coordinate lookup into promoter lookup"
		assert not self._coordlookup is None
		self._poverlaplookup = GenomicIntervalFinder("Gene Promoters")
		self._poverlaplookup.initPromotersFromCoordLookup(self._coordlookup, pbounds)

	def lookupName(self, gid):
		assert not self._coordlookup is None
		return self._coordlookup.getName(gid)

	def getGeneIDsByInterval(self, chrom, start, end):
		if self._overlaplookup is None:
			self.initOverlapsFromCoords()
		return self._overlaplookup.lookupInterval(chrom, start, end)

	def getPromotersByInterval(self, chrom, start, end, promoterbounds):
		assert self._poverlaplookup is not None
		return self._poverlaplookup.lookupInterval(chrom, start, end)

	def promoterOverlap(self, geneid, promoterbounds, elementgif):
		"""
		looks to see whether a geneid has any of the given elements in its promoter
		geneid -- a gene id (str)
		promoterbounds -- pair (tuple) of integers (upstream, downstream)
		elementgif -- a GenomicIntervalFinder object with the elements to lookup
		"""
		chrom, start, end = self._coordlookup.getCoord(geneid)
		strand = self._coordlookup.getStrand(geneid)
		if strand == "1":
			promoterstart = start - promoterbounds[0]
			promoterend = start + promoterbounds[1]
		if strand == "-1":
			promoterstart = start - promoterbounds[1]
			promoterend = start + promoterbounds[0]
		return elementgif.getOverlappingIntervalCount(chrom, promoterstart, promoterend)

if __name__ == "__main__":
	genefile = "/home/pzs/genebuilds/genelists/human-genes-ncbi36.csv"
	gfsettings = {	"id"	:	0,
					"strand":	2,
					"chrom"	:	3,
					"start"	:	4,
					"end"	:	5	}
	gl = GeneLookup()
	gl.initCoordsFromFile(genefile, gfsettings)
	gl.initOverlapsFromCoords()
	gl.initPromoterOverlapsFromCoords((500, 1000))
