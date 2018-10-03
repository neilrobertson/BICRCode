#!/usr/bin/env python

import csv
import os
import re

from cintervaltree import *
from bp import breakpoint

verbose = False

def isChrom(chrom):
	chrom = chrom.replace("chr", "")
	try:
		chromi = int(chrom)
		return chromi <= 23
	except ValueError:
		if chrom.lower() in [ "x", "y" ]:
			return True
		else:
			return False

wigsettings = {	"chrom"	:	0,
				"start"	:	1,
				"end"	:	2,
				"delimiter"	:	"\t"	}

class GenomicIntervalFinder:
	def __init__(self, name, padding=0):
		self.name = name
		# indexed by chromosome
		self._trees = {}
		# the value stored in each interval tree node, 
		# this can be overridden by specific values on initialisation
		self._defaultval = None
		# if we want to find hits within a tolerance, use the padding arg
		# this will increase the size of the interval by <padding> in both directions
		self._padding = padding

	#####
	# inits

	# this currently allows multiple calls, 
	# so you can fold several files into the Finder.
	# this behaviour is used by self.initFromGlob()
	# filesettings must have:
	#	delimiter
	# also EITHER
	#	coord
	# OR
	# 	chrom
	#	start
	#	end
	# OR
	#	chrom
	#	start
	#	padding
	# where padding is the distance from start where each interval should run
	# optional:
	#	value - a value at this interval
	#	end OR padding - end is the end of the interval. 
	#		padding means the end is start+padding
	#
	def initFromFile(self, filename, filesettings):
		print "initialising GenomicIntervalFinder", self.name, "from file", filename
		try:
			if "value" in filesettings:
				valueindex = filesettings["value"]
			if "padding" in filesettings:
				assert "end" not in filesettings
				padvalue = filesettings["padding"]
			if "end" in filesettings:
				assert "padding" not in filesettings
				endindex = filesettings["end"]
			if "coord" in filesettings:
				assert "start" not in filesettings
				assert "end" not in filesettings
				assert "chrom" not in filesettings
				assert "padding" not in filesettings
				coordindex = filesettings["coord"]
			if "chrom" in filesettings:
				chromindex = filesettings["chrom"]
				startindex = filesettings["start"]
			delimiter = filesettings["delimiter"]
		except KeyError:
			print "not all indices present"
			raise
		badlines = 0
		reader = csv.reader(open(filename), delimiter=delimiter)
		for row in reader:
			if len(row) == 1:
				continue
			if "value" in filesettings:
				value = row[valueindex]
			else:
				value = self._defaultval
			if "coord" in filesettings:
				coord = row[coordindex]
				mo = re.match("(\w+):(\d+)-(\d+)", coord)
				chrom = mo.group(1)
				start = int(mo.group(2))
				end = int(mo.group(3))
			else:
				chrom = row[chromindex]
				start = int(row[startindex])
				if "end" in filesettings:
					end = int(row[endindex])
				elif "padding" in filesettings:
					end = start + padvalue
				else:
					assert False
			if not isChrom(chrom):
				if verbose:
					print "invalid chrom on line", row
				badlines += 1
				continue
			self._addInterval(chrom, start, end, value)
		print "bad lines:", badlines
			
	def initFromWiggle(self, wigpath):
		self.initFromFile(wigpath, wigsettings)

	# convenience function for reading wiggle files
	def initWigglesFromGlob(self, wigdir):
		self.initFromGlob(wigdir + "/*.txt", wigsettings)

	def initFromGlob(self, pattern, filesettings):
		from glob import glob
		allfiles = glob(pattern)
		# if we sort these, the output makes more sense
		allfiles.sort()
		for ifile in allfiles:
			self.initFromFile(ifile, filesettings)

	def initFromCoordLookup(self, coordlookup):
		for id_ in coordlookup:
			chrom, start, end = coordlookup.getCoord(id_)
			self._addInterval(chrom, start, end, id_)

	def initPromotersFromCoordLookup(self, coordlookup, promoterbounds):
		up, down = promoterbounds
		pwriter = csv.writer(open("promoters.csv", "a"), delimiter="\t")
		for id_ in coordlookup:
			chrom, start, end = coordlookup.getCoord(id_)
			strand = coordlookup.getStrand(id_)
			if strand == "1":
				lower = start - up
				upper = start + down
			elif strand == "-1":
				lower = end - down
				upper = end + up
			else:
				assert False, strand
			self._addInterval(chrom, lower, upper, id_)
			pwriter.writerow( [chrom, lower, upper, id_ ])

	# end inits
	#####

	#####
	# traversals - all items in the tree

	def clearValues(self):
		from clearintervaltree import ClearIntervalTree
		setzero = ClearIntervalTree(0)
		self.applyTraversal(setzero)

	def applyTraversal(self, travfunc):
		"""
		applies a traversal over each tree
		travfunc -- a function (or class with __call__ implemented)
					to be applied to each node in each tree
		if travfunc implements a getFinalCount method, 
		this will be used to set a value for each tree.
		this is, for example, for counting the number of sites with hits in each tree
		"""
		travresults = {}
		for treename in self._trees:
			travresults[treename] = self._applyTraversalOneTree(treename, travfunc)
		return travresults

	def applyTraversalByTree(self, travfuncs):
		"""
		applies a traversal over each tree
		with a specific function for each tree
		travfuncs -- a dictionary mapping tree names to functions (or classes with __call__ implemented)
					to be applied to each node in each tree
					travfuncs.keys() must have the same items as self._trees.keys()
		"""
		assert set(travfuncs.keys()) == set(self._trees.keys())
		travresults = {}
		for treename in self._trees:
			thistravfunc = travfuncs[treename]
			self._applyTraversalOneTree(treename, thistravfunc)

	def betweenIntervalOverlap(self, gif):
		"""
		finds where there are hits from another set of genomic intervals
		in "fragments" - the gaps between the intervals we are currently storing 
		uses treetraversalobjects.OtherElementFragmentOverlap()
		gif -- Another GenomicIntervalFinder object that contains 
			the things we're looking for in between our intervals
		"""
		from treetraversalobjects import OtherElementFragmentOverlap
		assert set(gif.getTreeNames()) == set(self.getTreeNames())
		oneaway = 0
		twoaway = 0
		for treename in self.getTreeNames():
			# we need the main tree to find the "fragments" - stretches of genome
			# between the intervals
			thistree = self.getTree(treename)
			# we need the other tree to see when there are elements in these fragments
			othertree = gif.getTree(treename)
			fragmentcounter = OtherElementFragmentOverlap(thistree, othertree)
			thistree.traverse(fragmentcounter)
			oneaway += fragmentcounter.getOneAwayCount()
			twoaway += fragmentcounter.getTwoAwayCount()
		return oneaway, twoaway

	# end traversals
	#####

	#####
	# simple lookup functions
	def lookupInterval(self, chrom, start, end):
		"""
		gets the interval objects values that overlap the given parameters
		chrom -- name of chromosome. Should start with "chr"
		start -- genomic start
		end -- genomic end
		"""
		ivals = self._findIntervalOverlaps(chrom, start, end)
		return [ ival.value for ival in ivals ]
		
	def getOverlappingIntervalCount(self, chrom, start, end):
		ivals = self._findIntervalOverlaps(chrom, start, end)
		return len(ivals)
		
	def incrementOverlappingIntervals(self, chrom, start, end):
		"""
		increments all the overlapping intervals
		makes sure the intervals it finds have integer values
		returns True if it made any increments
		"""
		ivals = self._findIntervalOverlaps(chrom, start, end)
		try:
			assert all([ type(ival.value) == type(1) for ival in ivals ])
		except AssertionError:
			breakpoint()
		for ival in ivals:
			ival.value += 1
		return (len(ivals) > 0)

	def dumpOverlapsAsTrack(self, chrom, start, end, filepath, tracktype):
		"""
		finds the overlaps and attempts to produce a track (wiggle or bed file)
		chrom -- name of chromosome. Should start with "chr"
		start -- genomic start
		end -- genomic end
		filename -- where to dump track
		tracktype -- on of [ "bed", "wiggle" ]
		"""
		assert tracktype in [ "bed", "wiggle" ]
		print "writing %s file to %s" % (tracktype, filepath)
		position = "%s:%s-%s" % (chrom, start, end)
		fh = open(filepath, "w")
		print >> fh, "browser position %s\r" % (position)
		if tracktype == "bed":
			print >> fh, 'track name="%s" description="%s overlaps"\r' % (self.name, self.name)
		elif tracktype == "wiggle":
			print >> fh, 'track type=wiggle_0 name="%s" description="%s overlaps" visibility=full\r' % (self.name, self.name)
		intervals = self._findIntervalOverlaps(chrom, start, end)
		intervals.sort()
		writer = csv.writer(fh, delimiter="\t")
		if tracktype == "bed":
			rows = [ (chrom, i.start, i.end) for i in intervals ] 
		elif tracktype == "wiggle":
			rows = [ (chrom, i.start, i.end, i.value) for i in intervals ] 
		writer.writerows(rows)
		
	# end lookup functions
	#####


	#####
	# accessor/mutator
	def setDefaultValue(self, value):
		self._defaultval = value
	
	def getTreeNames(self):
		return self._trees.keys()

	def getTree(self, treename):
		return self._trees[treename]

	def setPadding(self, pad):
		self._padding = pad

	def getPadding(self):
		return self._padding

	# end accessor/mutator
	#####

	#####
	# private

	def _addInterval(self, chrom, start, end, value):
		pad = self._padding
		if not chrom.startswith("chr"):
			chrom = "chr" + chrom
		thistree = self._getTree(chrom)
		ival = Interval(start, end, value)
		thistree.insert_interval(ival)

	def _findIntervalOverlaps(self, chrom, start, end):
		try:
			assert end >= start
		except AssertionError:
			breakpoint()
		if not chrom.startswith("chr"):
			chrom = "chr" + chrom
		thistree = self._getTree(chrom)
		if thistree.isEmpty():
#			print "no tree for chromosome", chrom
#			breakpoint()
			return []
		else:
			pad = self._padding
			ivals = thistree.find(start-pad, end+pad)
#			ivals = thistree.find(start, end)
		return ivals

	def _getTree(self, chrom):
		tree = self._trees.setdefault(chrom, IntervalTree())
		return tree

	def _applyTraversalOneTree(self, treename, travfunc):
		summary = None
		tree = self._trees[treename]
		tree.traverse(travfunc)
		if hasattr(travfunc, "getFinalCount"):
			summary = travfunc.getFinalCount()
			travfunc.reset()
		return summary

	# end private
	#####

# example usage: load a vstep file and look up a few intervals
if __name__ == "__main__":
	vstepfile = os.path.expanduser("~/mount/publicdata/CD4-alignments/H3K4me3.bed.sort-binned.csv")
	vstepsettings = {	
						"chrom"	:	0,
						"start"	:	1,
						"padding"	:	200,
						"value"	:	2,
						"delimiter"	:	"\t",
					}
	gif = GenomicIntervalFinder("vstep-lookup")
	gif.initFromFile(vstepfile, vstepsettings)
	gif.dumpOverlapsAsTrack("chr1", 800000, 900000, "track.wig.txt", "wiggle")
	
