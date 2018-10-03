from bp import breakpoint

# function object for traversing interval trees
# this one counts up the number of intervals with a count >= minval in their "value" member
class IntervalCount:
	def __init__(self, minval):
		self._overcount = 0
		self._allcount = 0
		self._minval = minval
		
	def __call__(self, iv):
		thisval = iv.get_interval().value
		self._allcount += thisval
		if thisval >= self._minval:
			self._overcount += 1

	def reset(self):
		self._overcount = 0
		self._allcount = 0

	def getOverCount(self):
		return self._overcount
		
	def getAllCount(self):
		return self._allcount

	getFinalCount = getOverCount


# function object for traversing interval trees
# this one takes another set of interval trees for a different element (eg. CCACC)
# it counts the number of non-zero intervals in the main tree that have a CCACC 
# between it and neighbouring sites. The sections between sites are called "fragments"
class OtherElementFragmentOverlap:
	def __init__(self, basetree, elementtree):
		self._basetree = basetree
		self._elementtree = elementtree
		self._oneawaycount = 0
		self._twoawaycount = 0
		
	def _oneAwayHits(self, interval):
		try:
			leftboundary = self._basetree.before_interval(interval, max_dist=1e9)[0].start
			rightboundary = self._basetree.after_interval(interval, max_dist=1e9)[0].end
		except IndexError:
			breakpoint()
		return len(self._elementtree.find(leftboundary, rightboundary))

	def _twoAwayHits(self, interval):
		try:
			leftone = self._basetree.before_interval(interval, max_dist=1e9)[0]
			leftboundary = self._basetree.before_interval(leftone, max_dist=1e9)[0].start
			rightone = self._basetree.after_interval(interval, max_dist=1e9)[0]
			rightboundary = self._basetree.after_interval(rightone, max_dist=1e9)[0].end
		except IndexError:
			breakpoint()
		return len(self._elementtree.find(leftboundary, rightboundary))

	def __call__(self, iv):
		thisinterval = iv.get_interval()
		thisval = thisinterval.value
		if thisval > 0:
			self._oneawaycount += self._oneAwayHits(thisinterval)
			self._twoawaycount += self._twoAwayHits(thisinterval)
			
	def getOneAwayCount(self):
		return self._oneawaycount
		
	def getTwoAwayCount(self):
		return self._twoawaycount




# function object for traversing interval trees
# this one looks for intervals with value > 0
# and looks up the genes overlapping this site based on a list of tolerances
# (how far away the gene is allowed to be)
# it uses a GeneLookup() object to do this (self._gf)
class GeneCount:
	def __init__(self, chrom, tolerancelist):
		self._chrom = chrom.replace("chr", "")
		self._overlapcounters = {}
		self._tolerancelist = tolerancelist
		self._hitcount = 0
		self._misscount = 0
		self._misssites = []
		self._gf = None
		
	def updateGeneCount(self, geneid, genename, tolerance, numoverlaps):
		if geneid not in self._overlapcounters:
			newcounter = OneGeneOverlapCounter(geneid, genename, self._tolerancelist)
			self._overlapcounters[geneid] = newcounter
		overlapcounter = self._overlapcounters[geneid]
		overlapcounter.addSiteOverlaps(tolerance, numoverlaps)

	# iv is an Interval object from cintervaltree
	def processOneSite(self, iv):
		thisval = iv.get_interval().value
		if thisval > 0:
			didhit = False
			start, end = iv.start, iv.end
			for tolerance in self._tolerancelist:
				thisstart = start - tolerance
				thisend = end + tolerance
				geneoverlaps = self._gf.getGeneIDsByInterval(self._chrom, thisstart, thisend)
				if len(geneoverlaps) > 0:
					didhit = True
					for geneid in geneoverlaps:
						genename = self._gf.lookupName(geneid)
						self.updateGeneCount(geneid, genename, tolerance, thisval)
			if didhit:
				self._hitcount += 1
			else:
				self._misscount += 1
				self._misssites.append(iv)

	__call__ = processOneSite

	def getHitAndMissCount(self):
		return self._hitcount, self._misscount

	def setGeneFinder(self, gf):
		self._gf = gf

	def findPromoterElementOverlaps(self, promoterelements, promoterbounds):
		for geneid in self._overlapcounters:
			overlapcounter = self._overlapcounters[geneid]
			elcount = self._gf.promoterOverlap(geneid, promoterbounds, promoterelements)
			overlapcounter.setPromoterOverlaps(elcount)

	def dump(self, writer, tolerancelist, geneexp, usepromoterelements=False):
		for geneid in self._overlapcounters:
			overlapcounter = self._overlapcounters[geneid]
			expsummary = geneexp.getSummary(geneid)
			row = [ geneid, overlapcounter.name, self._chrom ] + expsummary + overlapcounter.getOverlaps(tolerancelist) 
			if usepromoterelements:
				row.append(overlapcounter.getPromoterOverlaps())
			writer.writerow(row)
			
	def dumpMisses(self, writer):
		for site in self._misssites:
			row = [ self._chrom, site.start, site.end, site.get_interval().value ]
			writer.writerow(row)

# companion class to store geneoverlap information for GeneCount()
class OneGeneOverlapCounter:
	def __init__(self, id_, name, tolerancelist):
		self.id = id_
		self.name = name
		self._overlapsbytolerance = dict([ (t, 0) for t in tolerancelist ])
		self._sitecount = dict([ (t, 0) for t in tolerancelist ])
		self._peoverlaps = 0
		
	def addSiteOverlaps(self, tolerance, overlaps):
		self._overlapsbytolerance[tolerance] += overlaps
		self._sitecount[tolerance] += 1
		try:
			assert self._sitecount[tolerance] <= self._overlapsbytolerance[tolerance]
		except:
			breakpoint()

	def getOverlaps(self, tolerancenames):
		# alternating site count and overlap count
		counts = []
		for tolerance in tolerancenames:
			counts.append(self._sitecount[tolerance])
			counts.append(self._overlapsbytolerance[tolerance])
		return counts

	def setPromoterOverlaps(self, promoteroverlaps):
		self._peoverlaps = promoteroverlaps
		
	def getPromoterOverlaps(self):
		return self._peoverlaps

