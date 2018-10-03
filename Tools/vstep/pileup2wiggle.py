# takes a pileup file, 
# which will have the same data value for a number of contiguous genome coordinate,
# and merge these contiguous points together.

# I've got carried away with generators here. I hope it makes sense

import sys
import csv

def readPilePositions(pileupfile):
	reader = csv.reader(open(pileupfile), delimiter="\t")
	for row in reader:
		yield row

# regions is an iterator and must be contiguous and sorted
def getCompressedRegions(regions):
	firstregion = regions.next()
	chrom, position, data = firstregion
	currentchrom = chrom
	currentstart = position
	lastposition = int(position)
	currentdata = data
	for region in regions:
		chrom, position, data = region
		position = int(position)
		# end of old region. Write and update new region start
		if chrom != currentchrom or data != currentdata or position > lastposition + 1:
			yield (currentchrom, currentstart, lastposition + 1, currentdata)
			currentchrom = chrom
			currentstart = position
			currentdata = data
		lastposition = position

def convertOnePileup(infile, outfile, name):
	print "writing compressed wiggle into file", outfile
	fh = open(outfile, "w")
	print >> fh, 'track type=wiggle_0 name="%s" visibility="full" description="%s"\r' % (name, name)
	writer = csv.writer(fh, delimiter="\t")
	pilepositionsgen = ( readPilePositions( infile ) )
	writer.writerows(( getCompressedRegions( pilepositionsgen ) ) )

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: pileup2wiggle.py <file> <name>"
		sys.exit(1)
	filepath = sys.argv[1]
	name = sys.argv[2]
	outpath = filepath + ".wig.txt"
	convertOnePileup(filepath, outpath, name)
