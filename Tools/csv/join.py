# there is a shell join, but it requires sorted files.

# this one works a bit like setcompare.py
# give it two files with an index
# it gives a combined file by joining on the given columns

# the column selected must contain each entry only once
# perhaps later we could relax this restriction

import sys
import csv

# maps the chosen column to the whole row
# silently remove anything with a multiple mapping
def getOneMapping(filename, column):
	fh = open(filename)
	multi = set()
	try:
		# if the read size is too small, it can't identify the delimiter properly
		dialect = csv.Sniffer().sniff(fh.read(5000), delimiters=" ,\t")
		fh.seek(0)
		reader = csv.reader(open(filename), dialect=dialect)
		colmap = {}
		for row in reader:
			id_ = row[column]
			if id_ in multi:
				continue
			if id_ in colmap:
				del colmap[id_]
				multi.add(id_)
			colmap[id_] = row
	except csv.Error:
		# couldn't find a dialect. Maybe this is only one column (not technically csv)
		fh.seek(0)
		assert column == 0
		colmap = {}
		for line in fh:
			line = line.rstrip()
			ids_ = line
			if id_ in multi:
				continue
			if id_ in colmap:
				del colmap[id_]
				multi.add(id_)
			colmap[id_] = line
	return colmap

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "usage: join.py <file1> <f1-column> <file2> <f2-column> "
		sys.exit(1)
	file1name = sys.argv[1]
	file1col = int(sys.argv[2])
	file2name = sys.argv[3]
	file2col = int(sys.argv[4])
	file1map = getOneMapping(file1name, file1col)
	file2map = getOneMapping(file2name, file2col)
	writer = csv.writer(sys.stdout, delimiter="\t")
	for id_ in file1map:
		if id_ not in file2map:
			continue
		file1row = file1map[id_]
		file2row = file2map[id_]
		file2row.remove(id_)
		writer.writerow(file1row + file2row)
	
