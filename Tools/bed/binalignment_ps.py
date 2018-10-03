#!/usr/bin/env python

import csv
import os

from bp import breakpoint

def getChromEnds(filename):
	print "getting chromends from file", filename
	chromends = {}
	reader = csv.reader(open(filename), delimiter="\t")
	# strip off header
	reader.next()
	for row in reader:
		chromname = row[0]
		if not chromname.startswith("chr"):
			chromname = "chr" + chromname
		chromends[chromname] = int(row[1])
	return chromends

# this is not quite a plane sweep because we don't dump all the values into a list first
class BinSolexaAlignmentsPlaneSweep:
	def __init__(self, chromends):
		self._chromends = chromends
		
	def dumpBin(self, writer, chrom, binstart, count):
		assert count >= 0
		if chrom is None:
			breakpoint()
		writer.writerow([ chrom, binstart, count ])

	# input files must be sorted by chromosome and position!
	def processOneFile(self, filename, binwidth):
		reader = csv.reader(open(filename), delimiter="\t")
		outfile = filename + "-binned.csv"
		if os.access(outfile, os.R_OK):
			print "already present", outfile
			return
		print "dumping binned alignments into file", outfile
		writer = csv.writer(open(outfile, "w"), delimiter="\t")
		lastchrom = None
		bintop = binwidth
		bincount = 0
		nextbincount = 0
		seenchrom = set()
		for row in reader:
			if not row:
				continue
			chrom, start, end = row
			start = int(start)
			end = int(end)
			if chrom not in self._chromends:
				print "unknown chromosome", chrom
				continue
			else:
				thischromend = self._chromends[chrom]
			assert start <= thischromend and end <= thischromend, "alignment outside chrom %d-%d / %d" % (start, end, thischromend)
			assert end > start
			# new chromosome
			if chrom != lastchrom:
				assert chrom not in seenchrom
				seenchrom.add(lastchrom)
				if bincount > 0:
					self.dumpBin(writer, lastchrom, (bintop - binwidth) + 1, bincount)
				bincount = 0
				bintop = binwidth
			# this alignment overlaps the current bin 
			if start < bintop:
				assert start >= bintop - binwidth
				bincount += 1
				if end > bintop + 1:
					nextbincount += 1
					assert end < bintop + binwidth
			# this alignment does not overlap
			else:
				# dump the last bin
				if bincount > 0:
					self.dumpBin(writer, chrom, (bintop - binwidth) + 1, bincount)
				# if the last bins neighbour also had alignments in it
				if nextbincount > 0:
					# shift up to the next bin
					bintop += binwidth
					if start < bintop:
						bincount = nextbincount + 1
						if end > bintop + 1:
							assert end < bintop + binwidth
							nextbincount = 1
						else:
							nextbincount = 0
						continue
					else:
						self.dumpBin(writer, chrom, (bintop - binwidth) + 1, nextbincount)
						nextbincount = 0
				# find the next bin where this alignment will overlap
				while start >= bintop:
					bintop += binwidth
				if bintop >= thischromend:
					print "bintop above chromend!", bintop, thischromend
				bincount = 1
				if end > bintop + 1:
					assert end < bintop + binwidth
					nextbincount += 1
			lastchrom = chrom
			
	def processManyFiles(self, basedir, pattern, binwidth):
		from glob import glob
		from os.path import basename
		allfiles = glob(basedir + "/*" + pattern)
		allnames = []
		for alignfile in allfiles:
			print "working at file", alignfile
			name = basename(alignfile).replace(pattern, "")
			allnames.append(name)
			self.processOneFile(alignfile, binwidth)
		return allnames

if __name__ == "__main__":
	import sys
	pattern = ".bed.sort"
	if len(sys.argv) != 4:
		print "usage: binalignment_ps.py <species> <basedir> <binwidth>"
		sys.exit(1)
	species = sys.argv[1]
	basedir = sys.argv[2]
	binwidth = int(sys.argv[3])

	chromendfile = os.path.expanduser("~/mount/publicdata/genebuilds/" + species + "/chromends.txt")
	chromends = getChromEnds(chromendfile)
	bsaps = BinSolexaAlignmentsPlaneSweep(chromends)
	allnames = bsaps.processManyFiles(basedir, pattern, binwidth)
	bw = csv.writer(open(basedir + "/basewidth.txt", "w"), delimiter=" ")
	rows = [ (name, binwidth) for name in allnames ]
	bw.writerows(rows)
