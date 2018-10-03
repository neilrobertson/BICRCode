import sys
import csv
from random import randint

class RandomBedGenerator:
	def __init__(self, chromendsfile):
		self.chromends = self.getChromEnds(chromendsfile)

	def getChromEnds(self, chromfile):
		reader = csv.reader(open(chromfile), delimiter="\t")
		# strip header
		reader.next()
		return dict([ (row[0], int(row[1])) for row in reader ])

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

	def readGenerator(self, readlen, numreads):
		chromgen = ( self.randomChromGenerator() )
		for i in range(numreads):
			chrom = chromgen.next()
			chromend = self.chromends[chrom]
			start = randint(0, chromend - readlen)
			end = start + readlen
			if not chrom.startswith("chr"):
				chrom = "chr" + chrom
			yield (chrom, start, end)

	def writeBed(self, bedfile, readlen, numreads):
		print "writing random bed into file", bedfile
		writer = csv.writer(open(bedfile, "w"), delimiter="\t")
		for read in self.readGenerator(readlen, numreads):
			writer.writerow(read)

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "usage: randombed.py <chromendfile> <readlen> <numreads> <outfile>"
		sys.exit(1)
	chromendfile = sys.argv[1]
	readlen = int(sys.argv[2])
	numreads = int(sys.argv[3])
	outfile = sys.argv[4]
	rbg = RandomBedGenerator(chromendfile)
	rbg.writeBed(outfile, readlen, numreads)
