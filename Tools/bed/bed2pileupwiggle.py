import csv
import sys
import os

from bp import breakpoint

def bedGenerator(bedfile):
	reader = csv.reader(open(bedfile), delimiter="\t")
	for row in reader:
		try:
			assert row[0].startswith("chr")
			yield(( row[0], int(row[1]), int(row[2]) ))
		except:
			print "bad row:", row
			raise

def sortBed(bedfile):
	sortcmd = "sort -k 1.4,1.5 -k 2,2n"
	sortedfile = bedfile + ".sort"
	cmd = "%s %s > %s" % (sortcmd, bedfile, sortedfile)
	print "sorting file with command", cmd
	os.system(cmd)
	return sortedfile

class RowBuffer(object):
	def __init__(self, writer):
		self._buffer = None
		self._writer = writer
		
	def writeRowBuffered(self, row):
		if row[1] == row[2] or row[3] == 0:
			return
		if self._buffer is not None:
			if self.isAdjacent(row):
				self.mergeRows(row)
				return
			else:
				self.writeBuffered()
		self._buffer = row
	
	def isAdjacent(self, row):
		# chromosomes the same. End of buffer matches start of new row. Values the same.
		return (	row[0] == self._buffer[0] and 
					row[1] == self._buffer[2] and
					row[3] == self._buffer[3])
	
	def mergeRows(self, row):
		self._buffer[2] = row[2]
	
	def writeRow(self, row):
		self._writer.writerow(row)

	def writeBuffered(self):
		self.writeRow(self._buffer)

def writeRow(writer, row):
#  	if row[0] == "chr1" and row[1] == 4775710:
#  		breakpoint()
	writer.writerow(row)

# for each new point we hit
	# deal with end points still on our list that are lower:
		# write region to from current start to the listed end point
		# decrement read count
		# set start of current region to end of last
	# set the current region start to the start of the point
	# add a 
def bedToPileup(bedgen, outfile):
	print "writing into file", outfile
	writer = csv.writer(open(outfile, "w"), delimiter="\t")
	
	currentstart = None
	currentdecrements = []
	lastchrom = None
	currentcount = 0
	stitch = False
	
	rb = RowBuffer(writer)
	
	for item in bedgen:
		chrom, start, end = item

		assert currentcount >= 0

		if chrom != lastchrom:
			while currentdecrements:
				regend = currentdecrements.pop(0)
				outrow = [ lastchrom, currentstart, regend + 1, currentcount ]
				rb.writeRowBuffered(outrow)
				currentcount -= 1
				currentstart = regend + 1
			currentstart = None
			assert currentcount == 0

		while currentdecrements and currentdecrements[0] < start:
			regend = currentdecrements.pop(0)
			outrow = [ chrom, currentstart, regend + 1, currentcount ]
			rb.writeRowBuffered(outrow)
			currentcount -= 1
			currentstart = regend + 1
			
		outrow = [ chrom, currentstart, start, currentcount ]
		rb.writeRowBuffered(outrow)
		currentcount += 1
		currentstart = start
		stitch = False
		currentdecrements.append(end)
		
		lastchrom = chrom
	while currentdecrements:
		regend = currentdecrements.pop(0)
		outrow = [ lastchrom, currentstart, regend + 1, currentcount ]
		rb.writeRowBuffered(outrow)
		currentcount -= 1
		currentstart = regend + 1
	rb.writeBuffered()

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: bed2pileupwiggle.py <bedfile> <outdir>"
		sys.exit(1)
	bedfile = sys.argv[1]
	outdir = sys.argv[2]

	sortedbed = sortBed(bedfile)
	print "working on file", sortedbed
	wigfile = "%s/%s.wig.txt" % (outdir, os.path.basename(sortedbed))
	bedToPileup(bedGenerator(sortedbed), wigfile)
