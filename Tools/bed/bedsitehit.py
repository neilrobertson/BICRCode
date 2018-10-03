import sys

# assume wiggle files
def buildSiteTree(sitedir):
	sitepattern = sitedir + "/*.txt"
	sitefiles = glob(sitepattern)
	git = GenomeIntervalTree()
	for sitefile in sitefiles:
		print "working at file", sitefile
		fh = open(sitefile)
		fh.next()
		fh.next()
		fh.next()
		reader = csv.reader(fh, delimiter="\t")
		for row in reader:
			chrom, start, end, data = row
			start, end = int(start), int(end)
			git.insertInterval(chrom, start, end, data)
	return git
	
def countSiteHits(bedfile, git):
	reader = csv.reader(open(bedfile), delimiter="\t")
	totalcount = 0
	hitcount = 0
	for row in reader:
		if len(row) == 0:
			print "bad row", row
			continue
		chrom, start, end = row
		hits = git.getOverlappingIntervals(chrom, start, end)
		if len(hits) > 0:
			hitcount += 1
		totalcout += 1
	print "%s: total: %d site hits %d" % (bedfile, totalcount, hitcount)

if __name__ == "__main__":
	from glob import glob
	sitedir = "/var/www/sreenivas/gtacwiggles/mouse"
	if len(sys.argv) != 2:
		print "usage bedsitehit.py <beddir>"
		sys.exit(1)
	git = buildSiteTree(sitedir)
	beddir = sys.argv[1]
	bedpattern = "%s/%s.bed" % (beddir)
	bedfiles = glob(bedpattern)
	for bedfile in bedfiles
