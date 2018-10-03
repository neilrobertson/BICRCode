import csv
import sys
import array
from random import randint

from datastructures.genomeintervaltree import *
from datastructures.clearintervaltree import *
from bp import breakpoint

######
# optimisation code to store chrom data with a smaller memory footprint
chrommap = {	'1'	:	"chr1", 
				'2'	:	"chr2", 
				'3'	:	"chr3", 
				'4'	:	"chr4", 
				'5'	:	"chr5", 
				'6'	:	"chr6", 
				'7'	:	"chr7", 
				'8'	:	"chr8", 
				'9'	:	"chr9", 
				'0'	:	"chr10", 
				'a'	:	"chr11", 
				'b'	:	"chr12", 
				'c'	:	"chr13", 
				'd'	:	"chr14", 
				'e'	:	"chr15", 
				'f'	:	"chr16", 
				'g'	:	"chr17", 
				'h'	:	"chr18", 
				'i'	:	"chr19", 
				'j'	:	"chr20", 
				'k'	:	"chr21", 
				'l'	:	"chr22", 
				'm'	:	"chrM",
				'x'	:	"chrX", 
				'y'	:	"chrY"
	}				

def reverseMap(cmap):
	rmap = {}
	for key in chrommap:
		value = chrommap[key]
		rmap[value] = key
	return rmap

rchrommap = reverseMap(chrommap)

rvg = None

def getRandomVstepGenerator(chromendfile, binsize):
	global rvg
	if rvg == None:
		print "building vstep generator object"
		rvg = RandomVstepGenerator(chromendfile, binsize)
	else:
		print "clearinging existing vstep generator object"
		cleartreetraverse = ClearIntervalTree(0)
		rvg._traverseAll(cleartreetraverse)
	return rvg

class GatherTraverse:
	def __init__(self):
		self.currentchrom = None
		self.chroms = array.array("c")
		self.coords = array.array("i")
		self.points = array.array("i")

	def setChrom(self, chromname):
		self.chrom = rchrommap[chromname]

	def __call__(self, ivalobj):
		ival = ivalobj.get_interval()
		if ival.value > 0:
			self.chroms.append(self.chrom)
			self.coords.append(ival.start)
			self.points.append(ival.value)

class DumpTraverse:
	def __init__(self, filename):
		print "dumping into file", filename
		self.writer = csv.writer(open(filename, "w"), delimiter="\t")
		self.chrom = None
		
	def setChrom(self, chromname):
		self.chrom = chromname

	def __call__(self, ivalobj):
		ival = ivalobj.get_interval()
		if ival.value > 0:
			row = [ self.chrom, ival.start, ival.end, ival.value ]
			self.writer.writerow(row)

class RandomVstepGenerator:
	def __init__(self, chromendsfile, binsize):
		self.chromends = self.getChromEnds(chromendsfile)
		self.binsize = binsize
		self.git = self.makeAllBinsIntervalTree()

	def getChromEnds(self, chromfile):
		reader = csv.reader(open(chromfile), delimiter="\t")
		return dict([ (row[0], int(row[1])) for row in reader ])

	def makeAllBinsIntervalTree(self):
		print "making vstep interval tree"
		git = GenomeIntervalTree()
		for chromname in self.chromends:
			print "working at chromosome", chromname
			chromend = self.chromends[chromname]
			index = 1
			while index < chromend:
				nextindex = index + self.binsize
				git.insertInterval(chromname, index, nextindex, 0)
				index = nextindex
		return git

	# this generates them in perpetuity
	def randomChromGenerator(self):
		chroms = self.chromends.keys()
		chroms.sort()
		cumulative = long(0)
		boundaries = []
		for chromname in chroms:
			chromlen = self.chromends[chromname]
			cumulative = cumulative + chromlen
			boundaries.append(cumulative)
		lastcoord = boundaries[-1]
		while 1:
			position = randint(0, lastcoord)
			yielded = False
			for i in range(len(boundaries)):
				thisboundary = boundaries[i]
				if position <= thisboundary:
					yield chroms[i]
					yielded = True
					break
			assert yielded, str(position) + str(boundaries)

	def recordRead(self, read):
		chrom, start, end = read
		intervals = self.git.getOverlappingIntervals(chrom, start, end)
		for interval in intervals:
			interval.value = interval.value + 1

	def readGenerator(self, readlen, numreads):
		chromgen = ( self.randomChromGenerator() )
		for i in range(numreads):
			chrom = chromgen.next()
			chromend = self.chromends[chrom]
			start = randint(0, chromend - readlen)
			end = start + readlen
			yield (chrom, start, end)

	# this piece of refactoring isn't quite finished
	# because _buildAndTraverse sets the chromosome for the traversal object
	# this isn't always necessary
	def _traverseAll(self, traverseobj):
		for chromname in self.git.values.keys():
			print "traversing intervals from chrom", chromname
			it = self.git.values[chromname]
			it.traverse(traverseobj)

	def _buildAndTraverse(self, readlen, numreads, traverseobj):
		readgen = ( self.readGenerator(readlen, numreads) )
		for read in readgen:
			self.recordRead(read)
		for chromname in self.git.values.keys():
			print "traversing intervals from chrom", chromname
			it = self.git.values[chromname]
			traverseobj.setChrom(chromname)
			it.traverse(traverseobj)

	def buildAndDumpVstep(self, readlen, numreads, outfile):
		dt = DumpTraverse(outfile)
		self._buildAndTraverse(readlen, numreads, dt)
	
	def buildAndGatherVstep(self, readlen, numreads):
		gt = GatherTraverse()
		self._buildAndTraverse(readlen, numreads, gt)
		return gt

if __name__ == "__main__":
	if len(sys.argv) != 6:
		print "usage: randomVstep.py <chromendfile> <readlen> <numreads> <binsize> <outfile>"
		sys.exit(1)
	chromendfile = sys.argv[1]
	readlen = int(sys.argv[2])
	numreads = int(sys.argv[3])
	binsize = int(sys.argv[4])
	outfile = sys.argv[5]
	rvg = RandomVstepGenerator(chromendfile, binsize)
	rvg.buildAndDumpVstep(readlen, numreads, outfile)
