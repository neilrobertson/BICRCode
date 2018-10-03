# take a vstep file
# find the median value (3rd column)
# divide all values by that
# spit out

import sys
import csv

# I haven't bothered to make this memory efficient
# so it may not work on huge files
def normByMedian(vstepfile):
	from numpy import median, mean
	reader = csv.reader(open(vstepfile), delimiter="\t")
	writer = csv.writer(sys.stdout, delimiter="\t")
	allrows = [ row for row in reader ]
	values = [ int(row[2]) for row in allrows ]
	avg = mean(values)
#	avg = median(values)
	writer.writerows([ (row[0], row[1], float(row[2]) / avg) for row in allrows ])

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "usage: normvstepbymedian.py <vstepfile>"
		sys.exit(1)
	vstepfile = sys.argv[1]
	normByMedian(vstepfile)
